
import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind, fisher_exact, norm
from scipy.io import savemat
from utils import build_obj, build_exp_methy_obj, add_to_block, format_expression_block_data, build_exp_singlegene_obj
from requests import get_methy_data, get_block_data, get_gene_data, get_sample_counts

from urllib.request import urlopen
import json
import itertools
import math
import statistical_methods
import UI_function
import data_functions

def comp_req(start_seq, end_seq, chromosome, gene_name, measurements=None):
    genomic_data = data_functions(start_seq, end_seq, chromosome, gene_name, measurements)
    stat_methods = statistical_methods()
    block_types = list (measurements.where(measurements['defaultChartType'] == 'block').dropna())
    gene_types = list (measurements.where(measurements['defaultChartType'] == 'scatterplot').dropna())
    methylation_df = measurements.where(measurements['defaultChartType'] == 'line').dropna()
    methylation_types = list (methylation_df.where(methylation_df['datasourceId'] == 'timp2014_probelevel_beta').dropna())
    methylation_diff_types = list (methylation_df.where(methylation_df['datasourceId'] != 'timp2014_probelevel_beta').dropna())
    per_gene_ttest = stat_methods.ttest_expression_per_gene(measurements, genomic_data.get_gene_data(),
                                               chromosome, start_seq, end_seq)

    yield per_gene_ttest
    # block overlap percentage
    block_overlap = stat_methods.block_overlap_percent(measurements, genomic_data.get_block_data(),
                                          start_seq, end_seq)
    yield block_overlap

    methy_raw_diff = get_methy_data(start_seq, end_seq, chromosome,measurements)

    methy_diff_corr_res = stat_methods.methy_diff_correlation(methy_raw_diff, measurements)

    yield methy_diff_corr_res

    methy_raw = get_methy_data(start_seq, end_seq, chromosome,measurements)
    methy_corr_res = stat_methods.methy_correlation(methy_raw, measurements)

        # loop through normal/tumor of each tissue type
    for data_source_one, data_source_two in itertools.combinations(
            measurements, 2):
        type1 = data_source_one["id"]
        type2 = data_source_two["id"]
        if type1.split("_")[0] != type2.split("_")[0]:
            continue

        correlation_coefficient = pearsonr(methy_raw[type1], methy_raw[
            type2])

        print (type1, type2)
        data_range = {
            'attr-one': [min(methy_raw[type1]), max(methy_raw[type1])],
            'attr-two': [min(methy_raw[type2]), max(methy_raw[type2])]
        }
        corr_obj = build_obj('correlation', 'methylation',
                             'methylation', True, data_source_one,
                             data_source_two,
                             correlation_coefficient[0],
                             correlation_coefficient[1],
                             ranges=data_range)
        methy_corr_res.append(corr_obj)
    methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
                            reverse=True)
    yield methy_corr_res


    corr_list = []
    # pvalue_list = []
    for data_source_one, data_source_two in itertools.combinations(
            measurements, 2):
        exp1 = data_source_one['id']
        exp2 = data_source_two['id']

        if exp1 not in measurements.columns or exp2 not in measurements.columns:
            continue

        col_one = measurements[exp1]
        col_two = measurements[exp2]

        correlation_coefficient = pearsonr(col_one, col_two)
        corr_obj = build_obj('correlation', 'expression', 'expression',
                             True, data_source_one,
                             data_source_two, correlation_coefficient[0],
                             correlation_coefficient[1])
        corr_list.append(corr_obj)

        t_value, p_value = ttest_ind(col_one, col_two,
                                     equal_var=False)
    corr_list = sorted(corr_list, key=lambda x: x['value'],
                       reverse=True)
    yield corr_list

    # gene expression and block independency test
    ttest_block_exp = stat_methods.ttest_block_expression(data_functions.get_gene_data(), data_functions.get_block_data(),
                                             gene_types, block_types)
    yield ttest_block_exp

        # correlation between methylation and gene expression
    # with the same tissue type
    corr_methy_gene = stat_methods.expression_methy_correlation(data_functions.get_gene_data(), gene_types, methylation_types,
                                                   methy_raw)

    yield corr_methy_gene

        # correlation between methylation difference and gene expression
    # difference
    corr_methy_gene = stat_methods.expression_methydiff_correlation(data_functions.get_gene_data(), gene_types, methylation_diff_types,
                                                       methy_raw_diff)

    yield corr_methy_gene
