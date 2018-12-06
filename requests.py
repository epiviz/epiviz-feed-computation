# handle requests to get data from data provider in EpiViz
#
# Created in Sep. 2017 by Zhe Cui
#

import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind

import urllib2
import json
import itertools
import math


# http request to get data
def get_url_data(data_source, measurements=None, chromosome=None,
                 start_seq=None, end_seq=None, metadata=None):
    # construct url
    # sql_url = 'http://epiviz-dev.cbcb.umd.edu/api/?requestId=10&version=4&' \
    #           'action=getData&datasource=' + data_source
    # sql_url = 'http://localhost:5000/?requestId=10&version=4&action=getData' \
    #           '&datasource=' + data_source
    sql_url = 'http://epiviz.umgear.org/api/?requestId=1&version=4&action=getData&datasourceGroup=hrp_2018&datasource=' + data_source    
    if measurements is not None:
        sql_url += '&measurement='
        if type(measurements) is list:
            sql_url += ','.join(measurements)
        else:
            sql_url += measurements

    if chromosome is not None:
        sql_url += '&seqName=' + str(chromosome)

    if start_seq is not None:
        sql_url += '&start=' + str(start_seq)
    if end_seq is not None:
        sql_url += '&end=' + str(end_seq)

    if metadata is not None:
        sql_url += '&metadata[]=' + str(metadata)

    # get data
    req = urllib2.Request(sql_url)
    response = urllib2.urlopen(req)
    a = json.loads(response.read())
    # check if it's an error message in the response, if it is, there's nothing coming back, we should skip the computation
    if a['error'] is not None and len(a['error']) > 1:
        return ""

    url_data = a['data']

    if url_data['rows']['useOffset']:
        relative_to_absolute(url_data['rows']['values'])
    # resolve the relative start and end position of the response
    # print url_data
    return url_data


def get_block_data(start, end, chromosome, block_measurements):
    # get block data
    block_data = dict()
    for block_measurement in block_measurements:
        # data source of block data: timp2014_breast_blocks, etc.
        # measurements: timp2014_breast_blocks, etc.
        # example: http://localhost:5000/?requestId=10&version=4
        # &action=getData&datasource=gene_methylation&measurement=cancer
        # &seqName=chr11&start=31164187&end=91164187
        # get data from url, the data is json format.
        url_data = get_url_data(data_source=block_measurement["datasourceId"],
                                chromosome=chromosome, start_seq=start,
                                end_seq=end)
                                
        # if there is no block data for this tissue type, just skip it
        if not url_data:
            continue

        # only keep start and end fields,  remove other data columns
        block_data[block_measurement["id"]] = pd.DataFrame(columns=["start",
                                                                    "end"])
        block_data[block_measurement["id"]]["start"] = url_data['rows'][
            'values']["start"]
        block_data[block_measurement["id"]]["end"] = url_data['rows'][
            'values']["end"]

    return block_data


def get_methy_data(start_seq, end_seq, chromosome, methylation_measurements):
    # get methylation data
    methy_raw = dict()
    for methylation_measurement in methylation_measurements:
        methy_data_source = methylation_measurement["datasourceId"]
        methy_id = methylation_measurement["id"]
        url_data = get_url_data(data_source=methy_data_source,
                                measurements=methy_id,
                                chromosome=chromosome,
                                start_seq=start_seq, end_seq=end_seq)

        # if there is no methy data for this tissue type, just skip it
        if not url_data:
            continue

        # extract expected format here
        curr_methy = dict()
        curr_methy['start'] = url_data['rows']['values']['start']
        curr_methy['end'] = url_data['rows']['values']['end']
        curr_methy['value'] = url_data['values']['values'][methy_id]
        methy_df = pd.DataFrame(curr_methy)
        methy_binned = methy_data_binning(start_seq, end_seq, 1000, methy_id, methy_df)

        methy_raw[methylation_measurement["name"]] = methy_binned['average']
        methy_raw['start'] = methy_binned['start']
        methy_raw['end'] = methy_binned['end']
        print "one methy type done!"

    return pd.DataFrame(methy_raw)


def get_gene_data(start_seq, end_seq, chromosome, gene_measurements):
    # get gene expressions data
    expression_data = dict()
    # gene expression has the same source
    gene_exp_data_source = gene_measurements[0]["datasourceId"]
    gene_exp_measurements = []
    for gene_measurement in gene_measurements:
        gene_exp_measurements.append(gene_measurement["id"])

    exp_url_data = get_url_data(data_source=gene_exp_data_source,
                                measurements=gene_exp_measurements,
                                chromosome=chromosome, start_seq=start_seq,
                                end_seq=end_seq, metadata="gene")

    if not exp_url_data:
        return pd.DataFrame(expression_data)

    # extract expected format here
    expression_data['start'] = exp_url_data['rows']['values']['start']
    expression_data['end'] = exp_url_data['rows']['values']['end']
    expression_data['gene'] = exp_url_data['rows']['values']['metadata']['gene']
    for tissue_type in gene_exp_measurements:
        # if there is no expression data for this tissue type, just skip it
        # if expression_data[tissue_type] is None:
        #     continue
        expression_data[tissue_type] = exp_url_data['values'][
            'values'][tissue_type]

    return pd.DataFrame(expression_data)


def get_sample_counts(measurements, start_seq, end_seq, chromosome):
    # construct url
    sql_url = 'http://54.157.53.251/api/?requestId=11&version=5&action=getValues&datasourceGroup=umd&datasource=gene_expression_barcode_subtype_count&metadata[]=gene'

    gene_exp_measurements = []
    for gene_measurement in measurements:
        gene_exp_measurements.append(gene_measurement["id"])

    sql_url += '&measurement='
    if type(gene_exp_measurements) is list:
        sql_url += ','.join(gene_exp_measurements)
    else:
        sql_url += gene_exp_measurements

    if chromosome is not None:
        sql_url += '&seqName=' + str(chromosome)

    if start_seq is not None:
        sql_url += '&start=' + str(start_seq)
    if end_seq is not None:
        sql_url += '&end=' + str(end_seq)

    # get data
    req = urllib2.Request(sql_url)
    response = urllib2.urlopen(req)
    a = json.loads(response.read())

    if a["error"]:
        return ""

    return a['data']['values']['values']


# convert relative values of start/end sequence number to absolute
def relative_to_absolute(raw_data):
    raw_data['start'] = np.cumsum(raw_data['start'])
    raw_data['end'] = np.cumsum(raw_data['end'])


def methy_data_binning(start_seq, end_seq, bins, methy_id, methy_data):
    bin_size = (end_seq - start_seq) / bins

    bin_number = 0
    results = dict()
    results['start'] = []
    results['end'] = []
    results['average'] = []
    start = start_seq
    end = start + bin_size - 1

    while bin_number < bins:
        bin_number += 1
        

        methy_filtered = methy_data[((start <= methy_data.start) & (
            methy_data.start <= end)) | ((start <= methy_data.end)
                                             & (methy_data.end <= end))]
        
        mean = methy_filtered.loc[:, "value"].mean()
        if not math.isnan(mean):
            results['average'].append(mean)
        else:
            results['average'].append(0)
        results['start'].append(start)
        results['end'].append(end)

        start += bin_size
        end += bin_size

    return results
        
