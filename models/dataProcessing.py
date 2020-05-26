import sys
import os
import ujson
import numpy as np
import pandas as pd
# import dask.dataframe as dd
import requests
import seaborn as sns
from sklearn.linear_model import LinearRegression
from epivizfileserver.measurements import WebServerMeasurement
# sys.path.insert(0,'/c/Users/Kyle/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/kych_linux/Documents/epiviz-feed-computation/src/epivizfeedcompute')
# print(sys.path)
from epivizfeedcompute.stats import BaseStats, TtestBlock, TtestExp, Correlation, CorrelationGeneSignal, OverlapBlock
import pickle

def split_n_buckets(dataframe, n):
    if len(dataframe) == 0:
        dataframe = [0 for i in range(n)]
    means = []
    dfs = np.array_split(dataframe, n,axis = 0) 
#     print(dfs)
    for i in range(n):
        means.append(np.mean(dfs[i]))
    return means
data = []
mid = "lung_cancer"
start=1
end=166151839
signal_data = pd.read_csv('timp2014_probelevel_beta_202001091200.csv')[["chr", "start", "end", "lung_normal", "lung_cancer"]]
signal_data = signal_data.where((signal_data.start >= start) & (signal_data.end <= end))
signal_data = signal_data.where((signal_data.chr != 'chrY'))
signal_data = signal_data.dropna()
cpg_islands = pd.read_csv('cpg_islands_202001141145.csv')
cpg_islands = cpg_islands.where((cpg_islands.chr != 'chrY')).dropna()

dataF = pd.DataFrame() 
for i in range(1,23):
    query = 'http://54.157.53.251/api/?requestId=11&version=5&action=getValues&datasourceGroup=umd&datasource=gene_expression_barcode_subtype&measurement=lung___tumor,lung___normal&seqName=chr{}&start=1&end=166151839&genome=&_=1573233064747'.format(i)
    result = requests.get(query)
    # res = umsgpack.unpackb(result.content)
    res = result.json()

    data = res['data']

    if data['rows']['useOffset']:
        data['rows']['values']['start'] = np.cumsum(data['rows']['values']['start'])
        data['rows']['values']['end'] = np.cumsum(data['rows']['values']['end'])

    # convert json to dataframe
    records = {}

    for key in data['rows']['values'].keys():
        if key not in ["id", "strand", "metadata"]:
            records[key] = data['rows']['values'][key]

    for key in data['rows']['values']['metadata'].keys():
        records[key] = data['rows']['values']['metadata'][key]

    for key in data['values']['values'].keys():
        records[key] = data['values']['values'][key]
    temp = pd.DataFrame(records)
    dataF = pd.concat([dataF,temp])

dataF.head()
for i in ['X']:
    query = 'http://54.157.53.251/api/?requestId=11&version=5&action=getValues&datasourceGroup=umd&datasource=gene_expression_barcode_subtype&measurement=lung___tumor,lung___normal&seqName=chr{}&start=1&end=166151839&genome=&_=1573233064747'.format(i)
    result = requests.get(query)
    # res = umsgpack.unpackb(result.content)
    res = result.json()

    data = res['data']

    if data['rows']['useOffset']:
        data['rows']['values']['start'] = np.cumsum(data['rows']['values']['start'])
        data['rows']['values']['end'] = np.cumsum(data['rows']['values']['end'])

    # convert json to dataframe
    records = {}

    for key in data['rows']['values'].keys():
        if key not in ["id", "strand", "metadata"]:
            records[key] = data['rows']['values'][key]

    for key in data['rows']['values']['metadata'].keys():
        records[key] = data['rows']['values']['metadata'][key]

    for key in data['values']['values'].keys():
        records[key] = data['values']['values'][key]
    temp = pd.DataFrame(records)
    dataF = pd.concat([dataF,temp])

# run this version
# signal_data = dd.from_pandas(signal_data, npartitions=10)
# dataF = dd.from_pandas(data, npartitions=10)
in_gene = []
methy_upstream = []
methy_downstream = []
upstream = 1000
downstream = 3000

for index, row in dataF.iterrows():
#     print(signal_data.where((((row.start <= signal_data.start) & (signal_data.start <= row.end)) | ((row.start <= signal_data.end) & (signal_data.end <= row.end))) & (row.chr == signal_data.chr)).dropna()[mid])
    in_gene_buckets = split_n_buckets(signal_data.where((((row.start <= signal_data.start) & (signal_data.start <= row.end)) | ((row.start <= signal_data.end) & (signal_data.end <= row.end))) & (row.chr == signal_data.chr)).dropna()[mid].tolist(), 20)
    methy_upstream_buckets = split_n_buckets(signal_data.where((((row.start - upstream<= signal_data.start) & (signal_data.start <= row.end)) | ((row.start - upstream <= signal_data.end) & (signal_data.end <= row.end))) & (row.chr == signal_data.chr)).dropna()[mid].tolist(), 10)
    methy_downstream_buckets = split_n_buckets(signal_data.where((((row.start <= signal_data.start) & (signal_data.start <= row.end + downstream)) | ((row.start <= signal_data.end) & (signal_data.end <= row.end + downstream))) & (row.chr == signal_data.chr)).dropna()[mid].tolist(), 10)
    in_gene.append(in_gene_buckets)
    methy_upstream.append(methy_upstream_buckets)
    methy_downstream.append(methy_downstream_buckets)

print("done 1")
in_gene_islands = []
upstream_islands = []
downstream_islands = []

for index, row in dataF.iterrows():
    # gets number of rows
    cpg_island_count = cpg_islands.where((row.start <= cpg_islands.start) & (row.end >= cpg_islands.end) & (row.chr == cpg_islands.chr)).dropna().shape[0]
    methy_upstream_cpg_island_count = cpg_islands.where(((row.start - upstream<= cpg_islands.start) & (row.start >= cpg_islands.end) & (row.chr == cpg_islands.chr))).dropna().shape[0]
    methy_downstream_cpg_island_count = signal_data.where(((cpg_islands.start >= row.end)  & (cpg_islands.end <= row.end + downstream) & (row.chr == cpg_islands.chr))).dropna().shape[0]

    contains_cpg_island = False
    if cpg_island_count > 0:
        contains_cpg_island = True
    in_gene_islands.append((cpg_island_count, contains_cpg_island))
    upstream_islands.append(methy_upstream_cpg_island_count)
    downstream_islands.append(methy_downstream_cpg_island_count)
print("done 2")

pickle.dump( in_gene, open( "in_gene_lung_tumor_buckets_v2.p", "wb" ) )

pickle.dump( methy_upstream, open( "methy_upstream_lung_tumor_buckets_v2.p", "wb" ) )

pickle.dump( methy_downstream, open( "methy_downstream_lung_tumor_buckets_v2.p", "wb" ) )
pickle.dump( in_gene_islands, open( "in_gene_lung_cpg_islands.p", "wb" ) )
pickle.dump( upstream_islands, open( "methy_upstream_lung_cpg_islands.p", "wb" ) )
pickle.dump( downstream_islands, open( "methy_downstream_lung_cpg_islands.p", "wb" ) )

pickle.dump( dataF, open( "expression_lung_data.p", "wb" ) )