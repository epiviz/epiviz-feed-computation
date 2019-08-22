from flask import Flask, Response
from interface import computational_request
from epivizFeed import LondFDR
from flask_cache import Cache
from flask_sockets import Sockets
import time
import ujson
import logging
import os
import pandas as pd
import json
import numpy as np

app = Flask(__name__)
cache = Cache(app, config={'CACHE_TYPE': 'simple', 'CACHE_DEFAULT_TIMEOUT': 0})

sockets = Sockets(app)

@sockets.route("/getInfo")
def info(websocket):
    with open( os.getcwd() + "/epiviz.json", "r") as config_file:
        data = ujson.loads(config_file.read())
        computations = data["computations"]
        measurements = data["measurements"]
        info = data["dataSources"]

    mgroups = {}
    for m in measurements:
        if m["datasourceGroup"] not in mgroups.keys():
            mgroups[m["datasourceGroup"]] = []
        mgroups[m["datasourceGroup"]].append(m["name"])
        
    for m in mgroups.keys():
        mgroups[m] = np.unique(mgroups[m])

    websocket.send(ujson.dumps({"sources": info, "groups": mgroups, "computations": computations}))
    
@sockets.route('/getdata')
def feed(websocket):
    message = ujson.loads(websocket.receive())
    
    with open( os.getcwd() + "/epiviz.json", "r") as config_file:
        data = ujson.loads(config_file.read())
        computations = data["computations"]
        measurements = data["measurements"]
        pval_threshold = data["pval_threshold"]
    
    logging.info(message)
    
    data = message['data']
    start = data['start']
    end = data['end']
    chromosome = data['chr']
    gene_name = data['gene']
    seqID = message['seq']
    pval_disc = data['significant']
    pval_count = data['totalTests']
    logging.info ("parameters")

    key = chromosome + '-' + str(start) + '-' + str(end)
    cached = cache.get(key)
    pValFDR = LondFDR(R = pval_disc, N = pval_count)
    if cached:
        significantCount = 0
        for cr in cached:
            if cr[u"significance"]:
                significantCount += 1
        websocket.send(ujson.dumps(cached))
        websocket.send(ujson.dumps({"seq": seqID, "significant": significantCount, "totalTests": len(cached)}))
        return
    
    results = computational_request(start, end, chromosome, gene_name, 
                    measurements=measurements, computations=computations, pval_threshold=pval_threshold)
    cache_results = []
    logging.info (results)

 
    for result in results:
        logging.info ("send back!")
        logging.info (time.time())
        logging.info ("\n")
        # emit('returned_results', result)
        significant_result =pd.DataFrame()

        # print(result)


        if len(result) > 0:
            logging.info("batch is > 1, FDR")
            pvalues = result["pvalue"].tolist()
            cor_pval = pValFDR.batchLondStar(pval = pvalues)
            result["adjusted_testing_levels"] = cor_pval[1]
            result["significance"] = cor_pval[2]
            # significant_result = result.loc[result['significance'] == True]
            
            pval_sci = [np.format_float_scientific(x, precision=10) for x in pvalues]
            result["pvalue"] = pval_sci
            comp_res = result.to_json(orient='records')
            parse_res = json.loads(comp_res)
            cache_results.extend(parse_res)
            websocket.send(ujson.dumps(parse_res))

    cache.set(key, cache_results)
    websocket.send(ujson.dumps({"seq": seqID, "significant": pValFDR.R, "totalTests": pValFDR.N}))

if __name__ == "__main__":
    from gevent import pywsgi
    from geventwebsocket.handler import WebSocketHandler
    server = pywsgi.WSGIServer(('', 5001), app, handler_class=WebSocketHandler)
    logging.info("Server Starts!")
    server.serve_forever()
