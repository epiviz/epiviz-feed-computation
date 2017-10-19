# Helper functions for statistical computations
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


def add_to_list(block, expression, tissue_key, block_key, ind):
    key = tissue_key + '|' + block_key
    if key not in block:
        block[key] = []
    block[key].append(expression[tissue_key][ind])


def add_to_block(tissue_types, gene_block, gene_expression, block_type, ind):
    for tissue in tissue_types:
        add_to_list(gene_block, gene_expression,
                    tissue, block_type, ind)
        # if tissue != 'pancreas':
        #     add_to_list(gene_block, gene_expression, tissue +
        #                 '___normal', block_type, ind)

