from flask import Flask, Response
from computation_request import computation_request
from flask_cache import Cache
from flask_sockets import Sockets
import time

import ujson

app = Flask(__name__)
cache = Cache(app, config={'CACHE_TYPE': 'simple', 'CACHE_DEFAULT_TIMEOUT': 0})

sockets = Sockets(app)
# socketio = SocketIO(app)

# @socketio.on('get_data_event', namespace='/getdata')
@sockets.route('/getdata')
def feed(websocket):
    message = ujson.loads(websocket.receive())
    measurements = test_measurements()
    print message
    data = message['data']
    start = data['start']
    end = data['end']
    chromosome = data['chr']
    seqID = message['seq']
    print "parameters"
    key = chromosome + '-' + str(start) + '-' + str(end)
    cached = cache.get(key)
    if cached:
        websocket.send(ujson.dumps(cached))
        websocket.send(ujson.dumps(seqID))
        return
    results = computation_request(start, end,
                                  chromosome,
                                  measurements=measurements)
    cache_results = []
    print results
    for result in results:
        print "send back!"
        print time.time()
        print "\n"
        # emit('returned_results', result)
        cache_results.extend(result)
        websocket.send(ujson.dumps(result))

    cache.set(key, cache_results)

    websocket.send(ujson.dumps(seqID))


# just for testing purposes
def test_measurements(expression=True, block=True, methylation=True):
    measurements = []
    gene_types = ['breast___normal', 'breast___tumor', 'colon___normal',
                  'colon___tumor', 'lung___normal', 'lung___tumor',
                  'thyroid___normal', 'thyroid___tumor']

    tissue_types = ['breast', 'colon', 'thyroid', 'lung']
    if expression:
        for gene_type in gene_types:
            measurements.append({
                "id": gene_type,
                "name": gene_type,
                "type": "feature",
                "datasourceId": "gene_expression_barcode_subtype",
                "datasourceGroup": "gene_expression_barcode_subtype",
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "scatterplot",
                "annotation": None,
                "metadata": ["probe"]
            })

    if block:
        for tissue_type in tissue_types:
            measurements.append({
                "id": 'timp2014_' + tissue_type + '_blocks',
                "name": 'timp2014_' + tissue_type + '_blocks',
                "type": "feature",
                "datasourceId": 'timp2014_' + tissue_type + '_blocks',
                "datasourceGroup": 'timp2014_' + tissue_type + '_blocks',
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "block",
                "annotation": None,
                "metadata": ["probe"]
            })

    if methylation:
        for tissue_type in tissue_types:
            measurements.append({
                "id": tissue_type,
                "name": tissue_type,
                "type": "feature",
                "datasourceId": 'timp2014_collapsed_diff',
                "datasourceGroup": 'timp2014_collapsed_diff',
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "line",
                "annotation": None,
                "metadata": ["probe"]
            })

    return measurements


# road map measurements
def roadmap_measurements(expression=True, block=True, methylation=True):
    measurements = []
    gene_types = ['breast___normal', 'breast___tumor', 'colon___normal',
                  'colon___tumor', 'lung___normal', 'lung___tumor',
                  'thyroid___normal', 'thyroid___tumor']

    tissue_types = ['breast', 'colon', 'thyroid', 'lung']
    if expression:
        for gene_type in gene_types:
            measurements.append({
                "id": gene_type,
                "name": gene_type,
                "type": "feature",
                "datasourceId": "gene_expression_barcode_subtype",
                "datasourceGroup": "gene_expression_barcode_subtype",
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "scatterplot",
                "annotation": None,
                "metadata": ["probe"]
            })

    if block:
        for tissue_type in tissue_types:
            measurements.append({
                "id": 'timp2014_' + tissue_type + '_blocks',
                "name": 'timp2014_' + tissue_type + '_blocks',
                "type": "feature",
                "datasourceId": 'timp2014_' + tissue_type + '_blocks',
                "datasourceGroup": 'timp2014_' + tissue_type + '_blocks',
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "block",
                "annotation": None,
                "metadata": ["probe"]
            })

    if methylation:
        for tissue_type in tissue_types:
            measurements.append({
                "id": tissue_type,
                "name": tissue_type,
                "type": "feature",
                "datasourceId": 'timp2014_collapsed_diff',
                "datasourceGroup": 'timp2014_collapsed_diff',
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "line",
                "annotation": None,
                "metadata": ["probe"]
            })

    return measurements


if __name__ == "__main__":
    from gevent import pywsgi
    from geventwebsocket.handler import WebSocketHandler
    server = pywsgi.WSGIServer(('', 5001), app, handler_class=WebSocketHandler)
    print "Server Starts!"
    server.serve_forever()


