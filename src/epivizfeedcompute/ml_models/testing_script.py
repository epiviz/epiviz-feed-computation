import sys
import os
import ujson
import time
from epivizfileserver.measurements import WebServerMeasurement
sys.path.insert(0,'C:/Users/Kyle/Documents/Python Scripts/Envs/research_epiviz/efeedcomputation/src/epivizfeedcompute')

from EpivizModel import EpivizModel
import helper_functions

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
            tm = WebServerMeasurement(m['type'], m['id'], m['name'],  "http://54.157.53.251/api/", 
                    m['datasourceId'], m['datasourceGroup'], m['annotation'], m['metadata'])
            if tm.annotation is None:
                tm.annotation = {}
            tm.annotation["datatype"] = m["datatype"]
            config_measurements.append(tm)
            
#ESR1 chr6: 150204511 - 157531913
# chr1: 0 - 1170420
chrom = "chr1"
start = 1
end = 20000
# filter for overlap measuremnets
overlap_measurements = [m for m in config_measurements if m.mid in [ "lung_cancer", "lung___tumor", "cpg_islands"] ]
for m in overlap_measurements:
    if m.mid in ["lung___tumor", "lung___normal"]:
        m.datatype = 'expression'
    elif m.mid in ["cpg_islands"]:
        m.datatype = 'peak'
    else:
        m.datatype = 'signal'
    # print(m.get_data(chrom, start, end))

def test_model():
    # create instance of the class
    # create instance of the class
    # orig_stdout = sys.stdout
    e = end
    # f = open('predictBenchmarks.txt', 'w')
    # sys.stdout = f
    st = time.time()
    test = EpivizModel(overlap_measurements, '.\epiviz_model_1.model', helper_functions.get_data)
    result = test.predict(chrom,start,e,{})
    print("{}bp time:{}".format(e-start, time.time()-st))
    e*=10
    result = test.predict(chrom,start,e,{})
    print("{}bp time:{}".format(e-start, time.time()-start))
    e*=10
    result = test.predict(chrom,start,e,{})
    print("{}bp time:{}".format(e-start, time.time()-start))
    e*=10
    result = test.predict(chrom,start,e,{})
    print("{}bp time:{}".format(e-start, time.time()-start))
    # sys.stdout = orig_stdout
    # f.close() 
    # test = EpivizModel(overlap_measurements, '.\epiviz_model_1.model', helper_functions.get_data)
    # result = test.predict(chrom,start,end,{})
    # assert len(result)
test_model()