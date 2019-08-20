from flask import Flask, Response
#from old_feed.computation_request import computation_request
from interface import computational_request
# from comp_req import comp_req
from epivizFeed import LondFDR
from flask_cache import Cache
from flask_sockets import Sockets
import time
import ujson
import logging
import os
import pandas as pd
import json

app = Flask(__name__)
cache = Cache(app, config={'CACHE_TYPE': 'simple', 'CACHE_DEFAULT_TIMEOUT': 0})

sockets = Sockets(app)
# socketio = SocketIO(app)

# @socketio.on('get_data_event', namespace='/getdata')

@sockets.route('/getdata')
def feed(websocket):
    message = ujson.loads(websocket.receive())
    
    with open( os.getcwd() + "/epiviz.json", "r") as config_file:
        measurements = ujson.loads(config_file.read())
    
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
    
    results = computational_request(start, end, chromosome, gene_name, measurements=measurements)
    cache_results = []
    logging.info (results)
 
    for result in results:
        logging.info ("send back!")
        logging.info (time.time())
        logging.info ("\n")
        # emit('returned_results', result)
        significant_result =pd.DataFrame()

        if len(result) > 1:
            logging.info("batch is > 1, FDR")
            pvalues = result["pvalue"].tolist()
            cor_pval = pValFDR.batchLondStar(pval = pvalues)
            result["adjusted_testing_levels"] = cor_pval[1]
            result["significance"] = cor_pval[2]
            # significant_result = result.loc[result['significance'] == True]


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
    # data = ujson.dumps(test_measurements())
    # with open('epiviz.json', 'w') as outfile:
    #     json.dump(test_measurements(), outfile)
    logging.info("Server Starts!")
    server.serve_forever()
