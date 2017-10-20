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


# http request to get data
def get_url_data(data_source, measurements=None, chromosome=None,
                 start_seq=None, end_seq=None):
    # construct url
    sql_url = 'http://localhost:5000/?requestId=10&version=4&action=getData' \
              '&datasource=' + data_source
    if measurements is not None:
        sql_url += '&measurement='
        if type(measurements) is list:
            sql_url += ','.join(measurements)
        else:
            sql_url += measurements

    if chromosome is not None:
        sql_url += '&seqName=chr' + str(chromosome)

    if start_seq is not None:
        sql_url += '&start=' + str(start_seq)
    if end_seq is not None:
        sql_url += '&end=' + str(end_seq)

    # get data
    req = urllib2.Request(sql_url)
    response = urllib2.urlopen(req)
    a = json.loads(response.read())
    url_data = a['data']
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
        # extract expected format here
        methy_raw['start'] = url_data['rows']['values']['start']
        methy_raw['end'] = url_data['rows']['values']['end']

        methy_raw[methy_id] = url_data['values']['values'][methy_id]

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
                                end_seq=end_seq)
    # extract expected format here
    expression_data['start'] = exp_url_data['rows']['values']['start']
    expression_data['end'] = exp_url_data['rows']['values']['end']
    for tissue_type in gene_exp_measurements:
        expression_data[tissue_type] = exp_url_data['values'][
            'values'][tissue_type]

    return pd.DataFrame(expression_data)


# convert relative values of start/end sequence number to absolute
def relative_to_absolute(raw_data):
    raw_data['start'] = np.cumsum(raw_data['start'])
    raw_data['end'] = np.cumsum(raw_data['end'])
