from flask import Flask, Response
from computation_request import ttest_block_expression, computation_request
from requests import get_methy_data, get_block_data, get_gene_data
from comp_req import comp_req
from flask_cache import Cache
from flask_sockets import Sockets
from run.py import test_measurements
import ttest_block
import time
import ujson

measurements = test_measurements()
start = data['start']
end = data['end']
chromosome = data['chr']
gene_name = data['gene']
gene_types = []
block_types = []

methylation_types = []
methylation_diff_types = []

for measurement in measurements:
    data_obj = {
        "id": measurement["id"],
        "name": measurement["name"],
        "datasourceId": measurement["datasourceId"]
    }
    if measurement["defaultChartType"] == "scatterplot":
        gene_types.append(data_obj)
    elif measurement["defaultChartType"] == "block":
        block_types.append(data_obj)
    elif measurement["defaultChartType"] == "line":
        if measurement["datasourceId"] == "timp2014_probelevel_beta":
            methylation_types.append(data_obj)
        else:
            methylation_diff_types.append(data_obj)

print(ttest_block_expression(gene_types, get_gene_data(start, end, chromosome, gene_types), chromosome, start, end))
print("this")
test = ttest_block(gene_types, get_gene_data(start, end, chromosome, gene_types), chromosome, start, end)
print(test.ttest())

