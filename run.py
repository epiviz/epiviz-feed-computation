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
    measurements = chipseq_measurements()
    print message
    data = message['data']
    start = data['start']
    end = data['end']
    chromosome = data['chr']
    gene_name = data['gene']
    seqID = message['seq']
    print "parameters"
    key = chromosome + '-' + str(start) + '-' + str(end)
    cached = cache.get(key)
    if cached:
        websocket.send(ujson.dumps(cached))
        websocket.send(ujson.dumps(seqID))
        return
    results = computation_request(start, end, chromosome, gene_name,
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


# chipseq data
def chipseq_measurements(expression=True, block=True, methylation=True):
    measurements = []

    gene_ids = ["p1hc1", "p1hc2", "p1hc3", "p1sc3", "p1sc2", "p1sc1", "p6sc2", "p6sc3", "p6hc2", "p6hc1"]

    gene_names = ["P1Hc RNA sample 1", "P1Hc RNA sample 2", "P1Hc RNA sample 3", "P1Sc RNA sample 3", "P1Sc RNA sample 2", "P1Sc RNA sample 1", "P6Sc RNA sample 2", "P6Sc RNA sample 3", "P6Hc RNA sample 2", "P6Hc RNA sample 1"]

    block_ids = ["p1hc_h3k9ac_peaks_merged", "p1sc_h3k9ac_peaks_merged", "P1_Hc_fAtoh1_GFP_filt_narrowPeak", "P1_Sc_Lfng_GFP_ppr_narrowPeak","P6_Hc_fAtoh1_GFP_ppr_IDR0_narrowPeak", "P6_Sc_Lfng_GFP_ppr_narrowPeak"]
    block_names = ["P1Hc H3k27m3 ChipSeq Peaks", "P1Hc ATAC Peaks", "P1Sc ATAC Peaks", "P1Sc H3k27m3 ChipSeq Peaks", "P6Hc ATAC Peaks", "P6Sc ATAC Peaks", "P6Sc H3k27m3 ChipSeq Peaks","P1HC H3k9ac ChipSeq Peaks", "P1Sc H3k9ac ChipSeq Peaks", "P6Hc H3k9ac ChipSeq Peaks", "P6Sc H3k9ac ChipSeq Peaks"]
    # ATAC blocks name: "P1_Hc_fAtoh1_GFP_filt_narrowPeak", "P1_Sc_Lfng_GFP_ppr_narrowPeak",
    # "P6_Hc_fAtoh1_GFP_ppr_IDR0_narrowPeak", "P6_Sc_Lfng_GFP_ppr_narrowPeak",
    # NOTE: "p6hc_h3k9ac_peaks_merged", "p6sc_h3k9ac_peaks_merged" is not available in the database for now, will add that later

    aggregation_datasource_groups = ["P1_Hc_atac-signal", "P1_Sc_atac-signal", "P6_Sc_atac-signal", "P6_Hc_atac-signal", "P6_Hc_atac-signal", "P6_Sc_atac-signal", "P1_Sc_atac-signal", "P1_Hc_atac-signal"]
    aggregation_data_names = ["P6Hc log ATAC signal", "P6Sc log ATAC signal", "P1Sc log ATAC signal", "P1Hc log ATAC signal",]
    aggregation_ids = ["lscore", "lscore", "lscore", "lscore"]

    if expression:
        for gene_id, gene_name in zip(gene_ids, gene_names):
            measurements.append({
                "id": gene_id,
                "name": gene_name,
                "type": "feature",
                "datasourceId": "rna_expression_merge",
                "datasourceGroup": "rna_expression_merge",
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "scatterplot",
                "annotation": None,
                "metadata": ["gene", "entrez", "gene_id"]
            })

    if block:
        for block_id, block_name in zip(block_ids, block_names):
            measurements.append({
                "id": block_id,
                "name": block_name,
                "type": "range",
                "datasourceId": block_id,
                "datasourceGroup": block_id,
                "dataprovider": "umd",
                "formula": None,
                "defaultChartType": "block",
                "annotation": None,
                "metadata": []
            })

    if methylation:
        # for methylation difference
        for datasource_group, name, idx in zip(aggregation_datasource_groups, aggregation_data_names, aggregation_ids):
            measurements.append({
                "id": idx,
                "name": name,
                "type": "feature",
                "datasourceId": datasource_group,
                "datasourceGroup": datasource_group,
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
