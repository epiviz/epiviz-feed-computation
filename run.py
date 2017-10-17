from flask import Flask, render_template
from computation_request import computation_request
import json
app = Flask(__name__)


@app.route('/chr=<chromosome>&start=<start>&end=<end>', methods=['GET', 'OPTIONS'])
def feed(start, end, chromosome):
    gene_types = ['breast___normal', 'breast___tumor', 'colon___normal',
                  'colon___tumor', 'lung___normal', 'lung___tumor',
                  'thyroid___normal', 'thyroid___tumor']

    tissue_types = ['breast', 'colon', 'thyroid', 'lung']
    results = computation_request(start, end, chromosome,
                                  gene_types=gene_types,
                                  block_types=tissue_types,
                                  methylation_types=tissue_types)
    print 'finished!'
    return json.dumps(results)
    # return results
    # return render_template('feed.html', results=results)


@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    return response

if __name__ == '__main__':
    app.run(debug=True, host='127.0.0.1', port=5001)


