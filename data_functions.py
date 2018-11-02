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


class data_funcations:

	def __init__(self, start_seq, end_seq, chromosome, gene_name, measurements=None):
		self.start_seq = start_seq
		self.end_seq = end_seq
		self.chromosome = chromosome
		self.gene_name = gene_name
		self.measurements = measurements


	def gene_data():

    	return get_gene_data(self.start_seq, self.end_seq, self.chromosome, self.measurements)

    def block_data():
    			
    	return get_block_data(self.start_seq, self.end_seq, self.chromosome, self.measurements)

    def methylation_diff():	
		return get_methy_data(start_seq, end_seq, chromosome,methylation_diff_types)

    def methylation():
    	return get_methy_data(start_seq, end_seq, chromosome,methylation_diff_types)

   	def methylation_correlation(data_source_one,data_source_two):
   		methy_raw = get_methy_data(start_seq, end_seq, chromosome,
                                   methylation_types)
        methy_corr_res = methy_correlation(methy_raw, methylation_diff_types)

        # loop through normal/tumor of each tissue type
        for data_source_one, data_source_two in itertools.combinations(
        	methylation_diff_types, 2):
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
        return methy_corr_res
  	
  	def expression_correlation():
        corr_list = []
        # pvalue_list = []
        for data_source_one, data_source_two in itertools.combinations(
        	gene_types, 2):
        exp1 = data_source_one['id']
        exp2 = data_source_two['id']

        if exp1 not in expression_data.columns or exp2 not in expression_data.columns:
        	continue

        	col_one = expression_data[exp1]
        	col_two = expression_data[exp2]

        	correlation_coefficient = pearsonr(col_one, col_two)
        	corr_obj = build_obj('correlation', 'expression', 'expression',
        		True, data_source_one,
        		data_source_two, correlation_coefficient[0],
        		correlation_coefficient[1])
        	corr_list.append(corr_obj)

        	t_value, p_value = ttest_ind(col_one, col_two,
        		equal_var=False)
            # ttest_obj = build_obj('t-test', 'expression', 'expression', True,
            #                       data_source_one, data_source_two, t_value,
            #                       p_value)
            # pvalue_list.append(ttest_obj)

        # pvalue_list = sorted(pvalue_list, key=lambda x: x['value'],
        #                      reverse=True)
        # yield pvalue_list

        corr_list = sorted(corr_list, key=lambda x: x['value'],
        	reverse=True)
        return corr_list
