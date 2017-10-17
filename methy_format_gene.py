import numpy as np
import pandas as pd
import urllib2
import json


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
        sql_url += '&seqName=' + str(chromosome)

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


# convert relative values of start/end sequence number to absolute
def relative_to_absolute(raw_data):
    raw_data['start'] = np.cumsum(raw_data['start'])
    raw_data['end'] = np.cumsum(raw_data['end'])


def methy_format_gene(start_seq, end_seq, chromosome, measurements,
                      downstream=3000, upstream=1000):
    exp_data = dict()
    methy_data = dict()
    corr_res = []
    exp_measurements = []
    methy_measurements = []
    exp_datasource = 'gene_expression_barcode_subtype'
    methy_datasource = 'timp2014_collapsed_diff'
    # get data from data provider
    for measurement in measurements:
        if 'expression' in measurement['name'].lower():
            exp_measurements.append(measurement['id'])
        else:
            methy_measurements.append(measurement['id'])

    url_methy = get_url_data(methy_datasource, methy_measurements, chromosome,
                             start_seq, end_seq)

    url_exp = get_url_data(exp_datasource, exp_measurements, chromosome,
                           start_seq, end_seq)

    # extract expected format here
    methy_data['start'] = url_methy['rows']['values']['start']
    methy_data['end'] = url_methy['rows']['values']['end']

    for tissue_type in methy_measurements:
        methy_data[tissue_type] = url_methy['values']['values'][tissue_type]

    # use pandas to get the data
    methy_pd = pd.DataFrame(methy_data)

    # gene expression data
    exp_data['start'] = url_exp['rows']['values']['start']
    exp_data['end'] = url_exp['rows']['values']['end']
    for tissue_type in exp_measurements:
        exp_data[tissue_type] = url_exp['values'][
            'values'][tissue_type]

    # use pandas to get the data
    exp_pd = pd.DataFrame(exp_data)
    methy_mean = pd.DataFrame(columns=methy_measurements )

    for i in range(len(exp_pd)):
        exp_start = exp_pd.iloc[i]['start'] - downstream
        exp_end = exp_pd.iloc[i]['end'] + upstream
        methy_filtered = methy_pd[(exp_start <= methy_pd.start) | (
            methy_pd.start <= exp_end) | (exp_start <= methy_pd.end)
            | (methy_pd.end <= exp_end)]
        # methy_filtered[methy_measurements].mean(axis=1)
        methy_mean = methy_mean.append(methy_filtered[
                                           methy_measurements].mean(),
                          ignore_index=True)

    # added_exp = exp_pd.add(methy_mean)
    for methy_measurement in methy_measurements:
        exp_pd[methy_measurement] = methy_mean[methy_measurement]
    return exp_pd
