import sys
import os
import ujson
from epivizfileserver.measurements import WebServerMeasurement
sys.path.insert(0,'C:/Users/Kyle/Documents/Python Scripts/Envs/research_epiviz/efeedcomputation/src/epivizfeedcompute')

from stats import BaseStats, TtestBlock, TtestExp, Correlation, CorrelationGeneSignal, OverlapBlock

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
    if m.mid == "colon___normal":
        m.annotation = {"sample_count": 16}
    else:
        m.annotation = {"sample_count": 486}
        
def test_ttest_exp():
    # create instance of the class
    test = TtestExp.TtestExp(overlap_measurements, 0.05)
    result = test.compute(chrom, start, end, {"datatype": "expression", "annotation":None})
    print(result)
    assert len(result)
test_ttest_exp()