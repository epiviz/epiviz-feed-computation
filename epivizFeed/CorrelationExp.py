import pandas as pd
import numpy as np
import math
import json
import itertools
from scipy.stats.stats import pearsonr, ttest_ind
from old_feed.utils import build_obj
from old_feed.UI_functions import format_exp_methy_output
from epivizFeed.StatMethod import StatMethod
from old_feed.data_functions import Gene_data, Methylation


class CorrelationExp(StatMethod):

    def __init__(self, measurements):
        super(CorrelationExp, self). __init__(measurements)
        self.gene_types = super(CorrelationExp, self).get_measurements_self("gene")

    def partion(self, part_type, group_one, group_two=None):
        #attributes other than tissues must be specifed in form [attr_name, (attr_val 1, attr_val 2)]
        # where attributes values must be distinct strings
        pair = None
        m = pd.DataFrame()
        for measurement in self.measurements:
            m = m.append(measurement, ignore_index=True)
        exp_measurements = m[m['name'].str.contains('Expression')]
        if part_type == None:
            pair = (exp_measurements, None)
        elif group_two is None:
            g_one = exp_measurements[exp_measurements['name'].str.contains(group_one)]
            g_two = exp_measurements[exp_measurements['name'].str.contains(group_one) == False]
            pair = (g_one, g_two)
        else:
            g_one = exp_measurements[exp_measurements['annotation'].str.contains(group_one)]
            g_two = exp_measurements[exp_measurements['annotation'].str.contains(group_two)]
            pair = (g_one, g_two)
        return pair

    def to_list_of_dict(self, group):
        ret_val = []
        for ele in group:
            for i in self.gene_types:
                if i["id"] == ele:
                    ret_val.append(i)

        return ret_val

    def unpack_params(self, additional):
        if additional is not None:
            partition_type = additional["partition_type"]
            group_one = additional["group_one"]
            group_two = additional["group_two"]
            grouping = additional["grouping"]
        else:
            partition_type = "condition"
            group_one = "normal"
            group_two = "tumor"
            grouping = "all_pairs"
        return partition_type, group_one, group_two, grouping

    def compute(self, chromosome, start, end, additional=None):
        part_type, g_one, g_two, grouping = self.unpack_params(additional)

        expression_data = Gene_data(start, end, chromosome, measurements=self.gene_types)
        corr_list = []

        group_one, group_two = self.partion(part_type, g_one)
        exp_group_one = Gene_data(start, end, chromosome, measurements=group_one.to_dict('records'))

        group_one = [c for c in exp_group_one.columns if "_" in c]
        group_one = self.to_list_of_dict(group_one)

        if part_type is not None:
            exp_group_two = Gene_data(start, end, chromosome, measurements=group_two.to_dict('records'))
            group_two = [c for c in exp_group_two.columns if "_" in c]
            group_two = self.to_list_of_dict(group_two)
            group_pairs = [(x, y) for x in group_one for y in group_two]
        else:
            group_pairs = itertools.combinations(group_one, 2)

        # all combinations of gene expressions
        # TODO (Kyle?): simplify the above code
        group_pairs = itertools.combinations(group_one + group_two, 2)
        # pvalue_list = []
        #for data_source_one, data_source_two in itertools.combinations(self.gene_types, 2):
        for data_source_one, data_source_two in group_pairs:
            exp1 = data_source_one['id']
            exp2 = data_source_two['id']

            if exp1 in expression_data.columns and exp2 in expression_data.columns:

                col_one = expression_data[exp1]
                col_two = expression_data[exp2]

                correlation_coefficient = pearsonr(col_one, col_two)
                corr_obj = build_obj('correlation', 'expression', 'expression',
                                     True, data_source_one,
                                     data_source_two, correlation_coefficient[0],
                                     correlation_coefficient[1])
                corr_list.append(corr_obj)

                t_value, p_value = ttest_ind(col_one, col_two, equal_var=False)

        corr_list = sorted(corr_list, key=lambda x: x['value'], reverse=True)
        corr_res = pd.Series(corr_list)
        corr_res = corr_res.apply(pd.Series)

        parse_res = corr_res
        # corr_res = corr_res.to_json(orient='records')
        # parse_res = json.loads(corr_res)

        return parse_res
