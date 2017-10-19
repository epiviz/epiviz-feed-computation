# This file computes correlations between different gene expression tracks.
#
# Created in Sep. 2017 by Zhe Cui
#

import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind
from utils import build_obj, add_to_block
from requests import get_methy_data, get_block_data, get_gene_data

import urllib2
import json
import itertools


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


# gene expression correlations and t tests
def correlation_computation(exp_data, exp_type):
    corr_list = []
    pvalue_list = []
    correlation_coefficients = exp_data[exp_type].corr()
    for exp1, exp2 in itertools.combinations(exp_type, 2):
        col_one = exp_data[exp1]
        col_two = exp_data[exp2]
        correlation_coefficient = correlation_coefficients.loc[exp1][exp2]
        t_value, p_value = ttest_ind(col_one, col_two,
                                     equal_var=False)
        corr_obj = build_obj('correlation', 'expression', 'expression',
                             True, exp1,
                             exp2, correlation_coefficient)
        corr_list.append(corr_obj)
        corr_list = sorted(corr_list, key=lambda x: x['value'],
                           reverse=True)
        ttest_obj = build_obj('ttest', 'expression', 'expression', True,
                              exp1, exp2, p_value)
        pvalue_list.append(ttest_obj)
        pvalue_list = sorted(pvalue_list, key=lambda x: x['value'],
                             reverse=True)
    return corr_list, pvalue_list


def ttest_block_expression(exp_data, block_data, exp_types, datasource_types):
    ttest_res = []
    gene_expression_block = dict()
    gene_expression_nonblock = dict()
    block_types = [datasource_type["id"] for datasource_type in
                   datasource_types]
    tissue_types = [exp_type["id"] for exp_type in exp_types]
    for block_type in block_types:

        gene_start_seq = exp_data['start']
        gene_end_seq = exp_data['end']
        # loop through gene expression
        for gene_ind in range(0, len(gene_start_seq)):
            gene_start = gene_start_seq[gene_ind]
            gene_end = gene_end_seq[gene_ind]
            if gene_in_block(block_data, block_type, gene_start, gene_end):
                add_to_block(tissue_types, gene_expression_block,
                             exp_data,
                             block_type, gene_ind)
            else:
                add_to_block(tissue_types, gene_expression_nonblock,
                             exp_data,
                             block_type, gene_ind)

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


def block_overlap_percent(data_sources, block_data):
    block_overlap = []
    tissue_types = [d['id'] for d in data_sources]
    for tissue_type_one, tissue_type_two in itertools.combinations(
            tissue_types, 2):
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

#
# def methy_correlation(tissue_types, methy_data):
#     methy_corr_res = []
#     correlation_coefficients = methy_data[tissue_types].corr()
#     for tissue_type1, tissue_type2 in itertools.combinations(tissue_types, 2):
#         correlation_coefficient = correlation_coefficients.loc[tissue_type1][
#             tissue_type2]
#         data_range = {
#             'attr-one': [min(methy_data[tissue_type1]), max(methy_data[
#                                                                   tissue_type1])],
#             'attr-two': [min(methy_data[tissue_type2]), max(methy_data[
#                                                                   tissue_type2])]
#         }
#         corr_obj = build_obj('correlation', 'methylation',
#                              'methylation', True, tissue_type1, tissue_type2,
#                              correlation_coefficient[0], ranges=data_range)
#         methy_corr_res.append(corr_obj)
#         methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
#                           reverse=True)
#
#     return methy_corr_res


def expression_methy_correlation(exp_data, datasource_gene_types,
                                 datasource_methy_types,
                                 methy_data, downstream=3000, upstream=1000):
    print "expression_methy_correlation"
    methy_types = [datasource_type["id"] for datasource_type in
                   datasource_methy_types]
    gene_types = [gene_type["id"] for gene_type in datasource_gene_types]
    methy_mean = pd.DataFrame(columns=methy_types)
    corr_res = []

    for i in range(len(exp_data)):
        exp_start = exp_data.iloc[i]['start'] - downstream
        exp_end = exp_data.iloc[i]['end'] + upstream
        methy_filtered = methy_data[((exp_start <= methy_data.start) & (
            methy_data.start <= exp_end)) | ((exp_start <= methy_data.end)
                                  & (methy_data.end <= exp_end))]
        methy_mean = methy_mean.append(methy_filtered[methy_types].mean(),
                                       ignore_index=True)

    # for now do not compute expression difference
    for tissue_type in gene_types:
        # use the difference of types (normal and tumor) for the same tissue
        # expression_type1 = tissue_type + '___normal'
        # expression_type2 = tissue_type + '___tumor'
        # expression_diff = np.subtract(exp_data['expression']
        #                               [expression_type1],
        #                               exp_data['expression'][
        #                                   expression_type2])
        expression_diff = exp_data[tissue_type]
        for methy_type in methy_mean:

            print tissue_type, methy_type
            correlation_coefficient = pearsonr(methy_mean[methy_type],
                                               expression_diff)
            print correlation_coefficient[0]

            # format the data into list of json objects for plots
            data = format_exp_methy_output(expression_diff, methy_mean[
                methy_type], tissue_type, methy_type)

            data_range = {
                'attr-one': [min(expression_diff),
                             max(expression_diff)],
                'attr-two': [min(methy_mean[methy_type]),
                             max(methy_mean[methy_type])]
            }
            corr_obj = build_obj('correlation', 'expression',
                                 'methylation', True, tissue_type, methy_type,
                                 correlation_coefficient[0],
                                 data=data, ranges=data_range)
            corr_res.append(corr_obj)
            corr_res = sorted(corr_res, key=lambda x: x['value'],
                              reverse=True)

    return corr_res


def format_exp_methy_output(attr1, attr2, type1, type2):
    data = []
    for exp, methy in zip(attr1, attr2):
        point = dict()
        point['expression_' + type1] = exp
        point['methylation_' + type2] = methy
        data.append(point)

    return data
    # methy_diff = dict()
    # for tissue_type in tissue_types:
    #     methy_diff[tissue_type] = []
    # methy_diff_by_ind = dict()
    #
    # gene_start_seq = exp_data['start']
    # gene_end_seq = exp_data['end']
    # methy_start_seq = methy_data['start']
    # methy_end_seq = methy_data['end']
    # gene_ind = 0
    # methy_ind = 0
    #
    # while gene_ind < len(gene_start_seq):
    #     if methy_ind >= len(methy_start_seq):
    #         break
    #     methy_start = methy_start_seq[methy_ind]
    #     methy_end = methy_end_seq[methy_ind]
    #     gene_start = gene_start_seq[gene_ind] - 3000
    #     gene_end = gene_end_seq[gene_ind] + 1000
    #
    #     # methy is before gene expression
    #     if gene_start > methy_end:
    #         methy_ind += 1
    #     # methy is after gene expression
    #     elif gene_end < methy_start:
    #         gene_ind += 1
    #     # methy is in between
    #     elif gene_start <= methy_end <= gene_end or \
    #                             gene_start <= methy_start <= gene_end or methy_start < gene_start < \
    #             methy_end or methy_start < gene_end < methy_end:
    #         if gene_ind not in methy_diff_by_ind:
    #             methy_diff_by_ind[gene_ind] = []
    #         methy_diff_by_ind[gene_ind].append(methy_ind)
    #         methy_ind += 1
    #
    # gene_ind = 0
    # while gene_ind < len(gene_start_seq):
    #     if gene_ind not in methy_diff_by_ind:
    #         for tissue_type in tissue_types:
    #             methy_diff[tissue_type].append(0)
    #     else:
    #         for tissue_type in tissue_types:
    #             target_sublist = [methy_data[tissue_type][i] for i in
    #                               methy_diff_by_ind[gene_ind]]
    #             methy_diff[tissue_type].append(np.mean(target_sublist))
    #     gene_ind += 1
    #
    # for tissue_type in tissue_types:
    #     for key, methy_value in methy_diff.items():
    #         # use the difference of types (normal and tumor) for the same tissue
    #         expression_type1 = tissue_type + '___normal'
    #         expression_type2 = tissue_type + '___tumor'
    #         expression_diff = np.subtract(exp_data['expression']
    #                                       [expression_type1],
    #                                       exp_data['expression'][
    #                                           expression_type2])
    #
    #         print tissue_type, key
    #         correlation_coefficient = pearsonr(methy_value,
    #                                            expression_diff)
    #         print correlation_coefficient[0]
    #
    #         # format the data into list of json objects for plots
    #         data = []
    #         for exp, methy in zip(expression_diff, methy_value):
    #             point = {}
    #             point['expression_' + tissue_type] = exp
    #             point['methylation_' + key] = methy
    #             data.append(point)
    #
    #         data_range = {
    #             'attr-one': [min(expression_diff),
    #                            max(expression_diff)],
    #             'attr-two': [min(methy_value), max(methy_value)]
    #         }
    #         corr_obj = build_obj('correlation', 'expression',
    #                              'methylation', True, tissue_type, key,
    #                              correlation_coefficient[0],
    #                              data=data, ranges=data_range)
    #         corr_res.append(corr_obj)
    #         corr_res = sorted(corr_res, key=lambda x: x['value'],
    #                           reverse=True)
    #         # corr_res.append(['expression ' + tissue_type, 'methylation '
    #         #                  'diff ' + key, correlation_coefficient[0]])
    #
    # return corr_res


def computation_request(start_seq, end_seq, chromosome, measurements=None):
    # skip pancreas for now, since we do not have that data
    # tissue_types = ['breast', 'colon', 'lung', 'thyroid']
    # start_seq = 3947953
    # end_seq = 7164991

    # extract data from measurements
    gene_types = []
    block_types = []
    methylation_types = []
    for measurement in measurements:
        data_obj = {
            "id": measurement["id"],
            "datasourceId": measurement["datasourceId"]
        }
        if measurement["defaultChartType"] == "scatterplot":
            gene_types.append(data_obj)
        elif measurement["defaultChartType"] == "block":
            block_types.append(data_obj)
        elif measurement["defaultChartType"] == "line":
            methylation_types.append(data_obj)

    return_results = []
    block_data = None
    methy_raw = None
    expression_data = None
    has_block = len(block_types) > 0
    has_methy = len(methylation_types) > 0
    has_gene = len(gene_types) > 0

    if has_block:
        block_data = get_block_data(start_seq, end_seq, chromosome, block_types)

        # block overlap percentage
        block_overlap = block_overlap_percent(block_types, block_data)
        return_results.extend(block_overlap)

    if has_methy:
        methy_types = [methylation_type["id"] for methylation_type in methylation_types]
        methy_raw = get_methy_data(start_seq, end_seq, chromosome,
                                   methylation_types)
        methy_corr_res = []
        # loop through every possible combinations of methylation
        for type1, type2 in itertools.combinations(methy_types, 2):
            correlation_coefficient = pearsonr(methy_raw[type1], methy_raw[
                type2])
            data_range = {
                'attr-one': [min(methy_raw[type1]), max(methy_raw[type1])],
                'attr-two': [min(methy_raw[type2]), max(methy_raw[type2])]
            }
            corr_obj = build_obj('correlation', 'methylation',
                                 'methylation', True, type1,
                                 type2,
                                 correlation_coefficient[0], ranges=data_range)
            methy_corr_res.append(corr_obj)
            methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
                                    reverse=True)

        return_results.extend(methy_corr_res)

    if has_gene:
        exp_types = [gene_type["id"] for gene_type in gene_types]
        expression_data = get_gene_data(start_seq, end_seq, chromosome,
                                        gene_types)

        # corr_exp, ttest_exp = correlation_computation(expression_data,
        #                                               gene_types)

        corr_list = []
        pvalue_list = []
        for exp1, exp2 in itertools.combinations(exp_types, 2):
            col_one = expression_data[exp1]
            col_two = expression_data[exp2]

            correlation_coefficient = pearsonr(col_one, col_two)
            corr_obj = build_obj('correlation', 'expression', 'expression',
                                 True, exp1,
                                 exp2, correlation_coefficient)
            corr_list.append(corr_obj)

            t_value, p_value = ttest_ind(col_one, col_two,
                                         equal_var=False)
            ttest_obj = build_obj('ttest', 'expression', 'expression', True,
                                  exp1, exp2, p_value)
            pvalue_list.append(ttest_obj)

        pvalue_list = sorted(pvalue_list, key=lambda x: x['value'],
                             reverse=True)
        corr_list = sorted(corr_list, key=lambda x: x['value'],
                           reverse=True)

        return_results.extend(corr_list)
        return_results.extend(pvalue_list)

    if has_gene and has_block:
        # gene expression and block independency test
        ttest_block_exp = ttest_block_expression(expression_data, block_data,
                                                 gene_types, block_types)
        return_results.extend(ttest_block_exp)

    if has_gene and has_methy:
        # correlation between methylation difference and gene expression
        # difference
        corr_methy_gene = expression_methy_correlation\
                         (expression_data, gene_types, methylation_types,
                          methy_raw)

        return_results.extend(corr_methy_gene)

    return return_results
