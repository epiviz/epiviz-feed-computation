from flask import Flask, stream_with_context, Response
from computation_request import computation_request
import json
app = Flask(__name__)


@app.route('/chr=<chromosome>&start=<start>&end=<end>', methods=['GET', 'OPTIONS'])
def feed(start, end, chromosome):
    # for testing purpose only

    measurements = test_measurements()
    # gene_types = ['breast___normal', 'breast___tumor', 'colon___normal',
    #               'colon___tumor', 'lung___normal', 'lung___tumor',
    #               'thyroid___normal', 'thyroid___tumor']
    #
    # tissue_types = ['breast', 'colon', 'thyroid', 'lung']
    results = computation_request(start, end, chromosome,
                                  measurements=measurements)

    # def generator():
    #     computation_request(start, end,
    #                         chromosome,
    #                         measurements=measurements)

    print 'finished!'
    # return stream_with_context(computation_request(start, end,
    #                         chromosome,
    #                         measurements=measurements))
    return json.dumps(results)
    return results
    # return render_template('feed.html', results=results)


@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    return response


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


if __name__ == '__main__':
    app.run(debug=True, host='127.0.0.1', port=5001)


