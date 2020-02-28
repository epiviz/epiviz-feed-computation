from .interface import computational_request
from .stats import LondFDR
from epivizfileserver.measurements import WebServerMeasurement
from epivizfileserver.client import EpivizClient

from flask import Flask, Response
from flask_caching import Cache
from flask_sockets import Sockets
from gevent import pywsgi
from geventwebsocket.handler import WebSocketHandler
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

def setup_app(server, file = None):
    global app
    app.config_file = file
    app.server = server
    app.measurements = []

    with open(app.config_file, "r") as config_file:
        data = ujson.loads(config_file.read())
        app.computations = data["computations"]
        app.info = data["dataSources"]
        app.pval_threshold = data["pval_threshold"]

        if data["measurements"] is not None:
            measurements = data["measurements"]

            for m in measurements: 
                tm = WebServerMeasurement(m['type'], m['id'], m['name'], app.server, 
                    m['datasourceId'], m['datasourceGroup'], m['annotation'], m['metadata'])
                if tm.annotation is None:
                    tm.annotation = {}
                tm.annotation["datatype"] = m["datatype"]
                app.measurements.append(tm)
        else:
            ec = EpivizClient(server)
            m = ec.get_measurements()
            app.measurements = m

    return app

def start_app(port=5001):
    global app
    server = pywsgi.WSGIServer(('', port), app, handler_class=WebSocketHandler)
    # print(app.config_file)
    logging.info("Server Starts!")
    server.serve_forever()

@sockets.route("/getInfo")
def info(websocket):
    with open(app.config_file, "r") as config_file:
        data = ujson.loads(config_file.read())
        computations = data["computations"]
        measurements = data["measurements"]
        pval_threshold = data["pval_threshold"]

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
    
    message = ujson.loads(websocket.receive())
    # with open(app.config_file, "r") as config_file:
    #     data = ujson.loads(config_file.read())
    #     computations = data["computations"]
    #     measurements = data["measurements"]
    #     pval_threshold = data["pval_threshold"]
    
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
    
    results = computational_request(chromosome, start, end, gene_name, 
                    measurements=app.measurements, computations=app.computations, pval_threshold=app.pval_threshold)
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

# if __name__ == "__main__":
#     from gevent import pywsgi
#     from geventwebsocket.handler import WebSocketHandler
#     server = pywsgi.WSGIServer(('', 5001), app, handler_class=WebSocketHandler)
#     logging.info("Server Starts!")
#     server.serve_forever()
start_app()