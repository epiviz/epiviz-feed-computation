import sys
import os
import ujson
from epivizfileserver.measurements import WebServerMeasurement
sys.path.insert(0,'C:/Users/Kyle/Documents/Python Scripts/Envs/research_epiviz/efeedcomputation/src/epivizfeedcompute')

from stats import BaseStats, TtestBlock, TtestExp, Correlation, CorrelationGeneSignal, OverlapBlock
# __author__ = "Jayaram Kancherla"
# __copyright__ = "Jayaram Kancherla"
# __license__ = "mit"​
# overlap_measurements = [
#     {
#         "id": "breast___normal",
#         "name": "Expression breast_normal",
#         "type": "feature",
#         "datasourceId": "gene_expression_barcode_subtype",
#         "datasourceGroup": "gene_expression_barcode_subtype",
#         "dataprovider": "umd",
#         "formula": None,
#         "datatype": "expression",
#         "defaultChartType": "scatterplot",
#         "annotation": None,
#         "metadata": ["probe"]
#     },
#     {
#         "id": "breast___tumor",
#         "name": "Expression breast_tumor",
#         "type": "feature",
#         "datasourceId": "gene_expression_barcode_subtype",
#         "datasourceGroup": "gene_expression_barcode_subtype",
#         "dataprovider": "umd",
#         "formula": None,
#         "datatype": "expression",
#         "defaultChartType": "scatterplot",
#         "annotation": None,
#         "metadata": ["probe"]
#     }
# ]
config_file = os.getcwd() + "/config.json"
config_measurements = []
with open(config_file, "r") as config_file:
    data = ujson.loads(config_file.read())
    # computations = data["computations"]
    # info = data["dataSources"]
    # pval_threshold = data["pval_threshold"]
    if data["measurements"] is not None:
        measurements = data["measurements"]
        for m in measurements: 
            config_measurements.append(
                WebServerMeasurement(m['type'], m['id'], m['name'], "http://54.157.53.251/api/", 
                m['datasourceId'], m['datasourceGroup'], m['annotation'], m['metadata']
                )
            )
#ESR1 chr6: 150204511 - 157531913
chrom = "chr6"
start = 150204511
end = 157531913
# filter for overlap measuremnets
overlap_measurements = [m for m in config_measurements if m.mid in ["colon___normal", "colon___tumor"] ]
for m in overlap_measurements:
    m.datatype = 'expression'
def test_corr():
    # create instance of the class
    test = Correlation.Correlation(overlap_measurements, 0.05)
    result = test.compute(chrom, start, end, {"datatype": "expression", "annotation":None})
    print(result)
    assert len(result)
test_corr()