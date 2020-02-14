import pandas as pd
import numpy as np

def split_n_buckets(dataframe, n):
    if len(dataframe) == 0:
        dataframe = [0 for i in range(n)]
    means = []
    dfs = np.array_split(dataframe, n,axis = 0) 
    for i in range(n):
        means.append(np.mean(dfs[i]))
    return means

def get_data(measurements, chr, start, end, params=None, upstream=1000, downstream=3000):
    for m in measurements:
        if m.datatype == "signal":
            signal_data = m.get_data(chr, start, end)[0]
          
            mid = m.mid
        elif m.datatype == "peak":
            cpg_islands = m.get_data(chr, start, end)[0]
    in_gene = split_n_buckets(signal_data.where((((start <= signal_data.start) & (signal_data.start <= end)) | ((start <= signal_data.end) & (signal_data.end <= end))) & (chr == signal_data.chr)).dropna()[mid].tolist(), 20)
    upstream = split_n_buckets(signal_data.where((((start - upstream<= signal_data.start) & (signal_data.start <= end)) | ((start - upstream <= signal_data.end) & (signal_data.end <= end))) & (chr == signal_data.chr)).dropna()[mid].tolist(), 10)
    downstream = split_n_buckets(signal_data.where((((start <= signal_data.start) & (signal_data.start <= end + downstream)) | ((start <= signal_data.end) & (signal_data.end <= end + downstream))) & (chr == signal_data.chr)).dropna()[mid].tolist(), 10)
    cpg_island_count = cpg_islands.where((start <= cpg_islands.start) & (end >= cpg_islands.end) & (chr == cpg_islands.chr)).dropna().shape[0]
    
    contains_islands = False
    if cpg_island_count > 0:
        contains_islands = True
    in_gene_islands = [cpg_island_count, contains_islands]
    
    in_gene_names = []
    upstream_names = []
    downstream_names = []
    for i in range(20):
        in_gene_names.append("in_gene {}".format(i))
        if i < 10:
            upstream_names.append("upstream {}".format(i))
            downstream_names.append("downstream {}".format(i))

    in_gene_df = pd.DataFrame([in_gene], columns = in_gene_names)
    upstream_df = pd.DataFrame([upstream], columns = upstream_names)
    downstream_df = pd.DataFrame([downstream], columns = downstream_names)
    cpg_islands_in_gene = pd.DataFrame([in_gene_islands],columns = ["cpg_count", "contains_cpg_islands"])
    
    X = pd.concat([in_gene_df,upstream_df,downstream_df,cpg_islands_in_gene], axis=1)
    X = X.fillna(0)
    X = np.expand_dims(X, axis=2)

    return X
