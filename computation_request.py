# This file computes correlations between different gene expression tracks.
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


# convert relative values of start/end sequence number to absolute
def relative_to_absolute(raw_data):
    raw_data['start'] = np.cumsum(raw_data['start'])
    raw_data['end'] = np.cumsum(raw_data['end'])


def gene_in_block(block_data, block_tissue, gene_start, gene_end):
    for block_start, block_end in zip(block_data[block_tissue]['start'],
                                      block_data[block_tissue]['end']):
        if block_start <= gene_end <= block_end:
            return True
        if block_start <= gene_start <= block_end:
            return True
        if gene_start <= block_start <= gene_end:
            return True
        if gene_start <= block_end <= gene_end:
            return True

    return False


def add_to_list(block, expression, tissue_key, block_key, ind):
    key = tissue_key + '|' + block_key
    if key not in block:
        block[key] = []
    block[key].append(expression[tissue_key][ind])


def add_to_block(tissue_types, gene_block, gene_expression, block_type, ind):
    for tissue in tissue_types:
        add_to_list(gene_block, gene_expression,
                    tissue + '___tumor', block_type, ind)
        if tissue != 'pancreas':
            add_to_list(gene_block, gene_expression, tissue +
                        '___normal', block_type, ind)


def sort_and_format(data, pos, reverse, split_char, type):
    formatted_results = []
    sort_rows = sorted(data, key=lambda x: x[pos], reverse=reverse)
    for element in sort_rows:
        formatted_results.append(type + ' between ' + element[0] + ' and ' +
                                 element[1] + ' is ' + str(element[2]) + '.')
    return formatted_results


def get_source_id(data_type, attribute_type):
    if data_type == 'expression':
        return 'gene_expression_barcode_subtype'
    elif data_type == 'block':
        return 'timp2014_' + attribute_type + '_blocks'
    elif data_type == 'methylation':
        return 'timp2014_collapsed_diff'


def build_obj(comp_type, data_one, data_two, show_chart, attr_one, attr_two,
              value, data=None, ranges=None):

    source_id_one = get_source_id(data_one, attr_one)
    source_id_two = get_source_id(data_two, attr_two)
    attr_one_range = [0, 1]
    attr_two_range = [0, 1]
    if ranges is not None:
        attr_one_range = ranges['attr-one']
        attr_two_range = ranges['attr-two']
    data_source = [{
        "id": source_id_one if data_one else attr_one,
        "name": attr_one,
        "type": "feature",
        "datasourceId": source_id_one,
        "datasourceGroup": source_id_one,
        "dataprovider": "umd",
        "formula": None,
        "defaultChartType": None,
        "annotation": None,
        "minValue": attr_one_range[0],
        "maxValue": attr_one_range[1],
        "metadata": ["probe"]
    }, {
        "id": source_id_two if data_two else attr_two,
        "name": attr_two,
        "type": "feature",
        "datasourceId": source_id_two,
        "datasourceGroup": source_id_two,
        "dataprovider": "umd",
        "formula": None,
        "defaultChartType": None,
        "annotation": None,
        "minValue": attr_two_range[0],
        "maxValue": attr_two_range[1],
        "metadata": ["probe"]
    }]
    target_obj = {
        'computation-type': comp_type,
        'data-type-one': data_one,
        'data-type-two': data_two,
        'show-chart': show_chart,
        'attribute-one': attr_one,
        'attribute-two': attr_two,
        'value': value,
        'data': data,
        'data-source': data_source
    }
    return target_obj


def correlation_computation(exp_data, exp_type):
    corr_list = []
    pvalue_list = []
    for exp1, exp2 in itertools.combinations(exp_data['expression'], 2):
        col_one = exp_data['expression'][exp1]
        col_two = exp_data['expression'][exp2]
        correlation_coefficient = pearsonr(col_one, col_two)
        t_value, p_value = ttest_ind(col_one, col_two,
                                     equal_var=False)
        corr_obj = build_obj('correlation', 'expression', 'expression',
                             True, exp1,
                             exp2, correlation_coefficient[0])
        corr_list.append(corr_obj)
        corr_list = sorted(corr_list, key=lambda x: x['value'],
                           reverse=True)
        ttest_obj = build_obj('ttest', 'expression', 'expression', True,
                              exp1, exp2, p_value)
        pvalue_list.append(ttest_obj)
        pvalue_list = sorted(pvalue_list, key=lambda x: x['value'],
                             reverse=True)
    return corr_list, pvalue_list


def ttest_block_expression(exp_data, block_data, tissue_types):
    ttest_res = []
    # construct block indicator
    block_indicator = dict()
    gene_expression_block = dict()
    gene_expression_nonblock = dict()
    for tissue_type in tissue_types:
        block_indicator[tissue_type] = []

        gene_start_seq = exp_data['start']
        gene_end_seq = exp_data['end']
        # loop through gene expression
        for gene_ind in range(0, len(gene_start_seq)):
            gene_start = gene_start_seq[gene_ind]
            gene_end = gene_end_seq[gene_ind]
            if gene_in_block(block_data, tissue_type, gene_start, gene_end):
                block_indicator[tissue_type].append(1)
                add_to_block(tissue_types, gene_expression_block,
                             exp_data[
                                 'expression'],
                             tissue_type, gene_ind)
            else:
                block_indicator[tissue_type].append(0)
                add_to_block(tissue_types, gene_expression_nonblock,
                             exp_data['expression'],
                             tissue_type, gene_ind)

    # calculate t test between block and non-block gene expression of the same
    # tissue type
    for key, gene_block_exp in gene_expression_block.items():
        if key not in gene_expression_nonblock:
            break
        gene_nonblock_exp = gene_expression_nonblock[key]
        print key
        t_value, p_value = ttest_ind(gene_block_exp, gene_nonblock_exp,
                                     equal_var=False)
        print p_value
        gene_type = key.split('|')[0]
        block_type = key.split('|')[1]
        ttest_obj = build_obj('ttest', 'expression', 'block', False,
                              gene_type, block_type, p_value)
        # ttest_res.append(['expression ' + gene_type, 'block ' + block_type,
        #                   p_value])
        ttest_res.append(ttest_obj)
        ttest_res = sorted(ttest_res, key=lambda x: x['value'],
                           reverse=True)

    return ttest_res


def block_overlap_percent(tissue_types, block_data):
    block_overlap = []
    # construct block indicator
    # overlap_ratio = dict()
    # block_intersection = dict()
    for tissue_type_one, tissue_type_two in itertools.combinations(
            tissue_types, 2):
        # overlap_ratio[tissue_type_one + '|' + tissue_type_two] = 0
        # block_intersection[tissue_type_one + '|' + tissue_type_two] = 0
        # loop through blocks for each tissue type
        # for blo in range(0, len(block_data[tissue_type_one][
        #                                       'start'])):
        block_tissue_one = block_data[tissue_type_one]
        block_tissue_two = block_data[tissue_type_two]
        block_one_ind = 0
        block_two_ind = 0
        block_one_len = len(block_tissue_one['start'])
        block_two_len = len(block_tissue_two['start'])

        if block_one_len < 1 or block_two_len < 1:
            continue
        overlap_region = []
        block_one_region = []
        block_two_region = []

        # calculate regions for each of the block tissues separately
        # union regions should be the sum of these regions minus overlap region
        for start, end in zip(block_tissue_one['start'], block_tissue_one[
            'end']):
            block_one_region.append(end - start)

        for start, end in zip(block_tissue_two['start'], block_tissue_two[
            'end']):
            block_two_region.append(end - start)

        while block_one_ind < block_one_len and block_two_ind < block_two_len:
            tissue_one_start = block_tissue_one['start'][block_one_ind]
            tissue_two_start = block_tissue_two['start'][block_two_ind]
            tissue_one_end = block_tissue_one['end'][block_one_ind]
            tissue_two_end = block_tissue_two['end'][block_two_ind]

            # there is an overlap
            if tissue_one_start <= tissue_two_start < tissue_one_end or \
               tissue_one_start < tissue_two_end <= tissue_one_end or \
               tissue_two_start <= tissue_one_start < tissue_two_end or \
               tissue_two_start < tissue_one_end <= tissue_two_end:
                overlap_region.append(min(tissue_two_end, tissue_one_end) - max(
                                     tissue_one_start, tissue_two_start))
                # union_region = max(tissue_two_end, tissue_one_end) - min(
                #                      tissue_one_start, tissue_two_start)
                # generally we can do this for either block
                block_one_ind += 1
            # block tissue two is larger
            elif tissue_two_start >= tissue_one_end:
                block_one_ind += 1
            # block tissue one is larger
            elif tissue_one_start >= tissue_two_end:
                block_two_ind += 1

        overlap = sum(overlap_region)
        union = sum(block_one_region) + sum(block_two_region) - overlap
        overlap_obj = build_obj('overlap', 'block', 'block', False,
                                tissue_type_one, tissue_type_two,
                                overlap * 1.0 / union)
        block_overlap.append(overlap_obj)

    block_overlap = sorted(block_overlap, key=lambda x : x['value'],
                           reverse=True)
    print 'overlap done!'
    return block_overlap


def methy_correlation(tissue_types, methy_data):
    methy_corr_res = []
    for tissue_type1, tissue_type2 in itertools.combinations(tissue_types, 2):
        correlation_coefficient = pearsonr(methy_data[tissue_type1],
                                           methy_data[tissue_type2])
        data_range = {
            'attr-one': [min(methy_data[tissue_type1]), max(methy_data[
                                                                  tissue_type1])],
            'attr-two': [min(methy_data[tissue_type2]), max(methy_data[
                                                                  tissue_type2])]
        }
        corr_obj = build_obj('correlation', 'methylation',
                             'methylation', True, tissue_type1, tissue_type2,
                             correlation_coefficient[0], ranges=data_range)
        methy_corr_res.append(corr_obj)
        methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
                          reverse=True)

    return methy_corr_res


def expression_methy_correlation(exp_data, tissue_types, methy_data):
    corr_res = []

    methy_diff = dict()
    for tissue_type in tissue_types:
        methy_diff[tissue_type] = []
    methy_diff_by_ind = dict()

    gene_start_seq = exp_data['start']
    gene_end_seq = exp_data['end']
    methy_start_seq = methy_data['start']
    methy_end_seq = methy_data['end']
    gene_ind = 0
    methy_ind = 0

    while gene_ind < len(gene_start_seq):
        if methy_ind >= len(methy_start_seq):
            break
        methy_start = methy_start_seq[methy_ind]
        methy_end = methy_end_seq[methy_ind]
        gene_start = gene_start_seq[gene_ind] - 3000
        gene_end = gene_end_seq[gene_ind] + 1000

        # methy is before gene expression
        if gene_start > methy_end:
            methy_ind += 1
        # methy is after gene expression
        elif gene_end < methy_start:
            gene_ind += 1
        # methy is in between
        elif gene_start <= methy_end <= gene_end or \
                                gene_start <= methy_start <= gene_end or methy_start < gene_start < \
                methy_end or methy_start < gene_end < methy_end:
            if gene_ind not in methy_diff_by_ind:
                methy_diff_by_ind[gene_ind] = []
            methy_diff_by_ind[gene_ind].append(methy_ind)
            methy_ind += 1

    gene_ind = 0
    while gene_ind < len(gene_start_seq):
        if gene_ind not in methy_diff_by_ind:
            for tissue_type in tissue_types:
                methy_diff[tissue_type].append(0)
        else:
            for tissue_type in tissue_types:
                target_sublist = [methy_data[tissue_type][i] for i in
                                  methy_diff_by_ind[gene_ind]]
                methy_diff[tissue_type].append(np.mean(target_sublist))
        gene_ind += 1

    for tissue_type in tissue_types:
        for key, methy_value in methy_diff.items():
            # use the difference of types (normal and tumor) for the same tissue
            expression_type1 = tissue_type + '___normal'
            expression_type2 = tissue_type + '___tumor'
            expression_diff = np.subtract(exp_data['expression']
                                          [expression_type1],
                                          exp_data['expression'][
                                              expression_type2])

            print tissue_type, key
            correlation_coefficient = pearsonr(methy_value,
                                               expression_diff)
            print correlation_coefficient[0]

            # format the data into list of json objects for plots
            data = []
            for exp, methy in zip(expression_diff, methy_value):
                point = {}
                point['expression_' + tissue_type] = exp
                point['methylation_' + key] = methy
                data.append(point)

            data_range = {
                'attr-one': [min(expression_diff),
                               max(expression_diff)],
                'attr-two': [min(methy_value), max(methy_value)]
            }
            corr_obj = build_obj('correlation', 'expression',
                                 'methylation', True, tissue_type, key,
                                 correlation_coefficient[0],
                                 data=data, ranges=data_range)
            corr_res.append(corr_obj)
            corr_res = sorted(corr_res, key=lambda x: x['value'],
                              reverse=True)
            # corr_res.append(['expression ' + tissue_type, 'methylation '
            #                  'diff ' + key, correlation_coefficient[0]])

    return corr_res


def get_block_data(start, end, chromosome, block_types):
    # get block data
    block_data = dict()
    for tissue_type in block_types:
        # data source of block data: timp2014_breast_blocks, etc.
        # measurements: timp2014_breast_blocks, etc.
        # example: http://localhost:5000/?requestId=10&version=4
        # &action=getData&datasource=gene_methylation&measurement=cancer
        # &seqName=chr11&start=31164187&end=91164187
        block_data_source = 'timp2014_' + tissue_type + '_blocks'
        # get data from url, the data is json format.
        url_data = get_url_data(data_source=block_data_source,
                                chromosome=chromosome, start_seq=start,
                                end_seq=end)
        block_data[tissue_type] = url_data['rows']['values']

    return block_data


def get_methy_data(start_seq, end_seq, chromosome, methylation_types):
    # get methylation data
    methy_data_source = 'timp2014_collapsed_diff'
    methy_raw = dict()
    url_data = get_url_data(methy_data_source, methylation_types, chromosome,
                            start_seq, end_seq)
    # extract expected format here
    methy_raw['start'] = url_data['rows']['values']['start']
    methy_raw['end'] = url_data['rows']['values']['end']

    for tissue_type in methylation_types:
        methy_raw[tissue_type] = url_data['values']['values'][tissue_type]

    return methy_raw


def get_gene_data(start_seq, end_seq, chromosome, gene_types):
    # get gene expressions data
    expression_data = dict()
    gene_exp_data_source = 'gene_expression_barcode_subtype'
    gene_exp_measurements = gene_types

    exp_url_data = get_url_data(gene_exp_data_source, gene_exp_measurements,
                                chromosome, start_seq, end_seq)
    # extract expected format here
    expression_data['start'] = exp_url_data['rows']['values']['start']
    expression_data['end'] = exp_url_data['rows']['values']['end']
    expression_data['expression'] = dict()
    for tissue_type in gene_exp_measurements:
        expression_data['expression'][tissue_type] = exp_url_data['values'][
            'values'][tissue_type]

    return expression_data


def computation_request(start_seq, end_seq, chromosome, gene_types=None,
                        block_types=None, methylation_types=None):
    # skip pancreas for now, since we do not have that data
    # tissue_types = ['breast', 'colon', 'lung', 'thyroid']
    chr_seq = chromosome
    # start_seq = 3947953
    # end_seq = 7164991

    return_results = []
    block_data = None
    methy_raw = None
    expression_data = None

    if block_types is not None:
        block_data = get_block_data(start_seq, end_seq, chromosome, block_types)

        # block overlap percentage
        block_overlap = block_overlap_percent(block_types, block_data)
        return_results.extend(block_overlap)

    if methylation_types is not None:
        methy_raw = get_methy_data(start_seq, end_seq, chromosome,
                                   methylation_types)

        corr_within_methy = methy_correlation(methylation_types, methy_raw)
        return_results.extend(corr_within_methy)

    if gene_types is not None:
        expression_data = get_gene_data(start_seq, end_seq, chromosome,
                                        gene_types)

        # iterate through the names of each columns
        # correlation between different gene expressions
        corr_exp, ttest_exp = correlation_computation(expression_data,
                                                  gene_types)
        return_results.extend(corr_exp)
        return_results.extend(ttest_exp)

    if gene_types is not None and block_types is not None:
        # gene expression and block independency test
        ttest_block_exp = ttest_block_expression(expression_data, block_data,
                                                 block_types)
        return_results.extend(ttest_block_exp)

    if gene_types is not None and methylation_types is not None:
        # correlation between methylation difference and gene expression
        # difference
        corr_methy_gene = expression_methy_correlation \
                         (expression_data, methylation_types, methy_raw)

        return_results.extend(corr_methy_gene)

    return return_results
