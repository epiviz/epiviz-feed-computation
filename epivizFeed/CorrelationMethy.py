import pandas as pd
import numpy as np
import math
import json
import itertools
from scipy.stats.stats import pearsonr, ttest_ind
from old_feed.utils import build_obj
from old_feed.UI_functions import format_exp_methy_output
from epivizFeed.StatMethod import StatMethod
from old_feed.data_functions import Methylation, Methylation_diff


class CorrelationMethy(StatMethod):

    def __init__(self, measurements, methy_type, pval_threshold):
        super(CorrelationMethy, self). __init__(measurements, pval_threshold)
        self.methylation_types = super(CorrelationMethy, self).get_measurements_self(methy_type)
        self.methy_type = methy_type

    def get_methy_data(self, start, end, chromosome):
        if self.methy_type == "methy":
            methy_data = Methylation(start, end, chromosome, measurements=self.methylation_types)
        else:
            methy_data = Methylation_diff(start, end, chromosome, measurements=self.methylation_types)
        return methy_data

    def partion(self, type, group_one, group_two=None):
        #attributes other than tissues must be specifed in form [attr_name, (attr_val 1, attr_val 2)]
        # where attributes values must be distinct strings
        pair = None
        m = pd.DataFrame()
        for measurement in self.measurements:
            m = m.append(measurement, ignore_index=True)
        keywords_in_measurements = 'Probe' if self.methy_type == 'methy' else 'Methylation'
        methy_measurements = m[m['name'].str.contains(keywords_in_measurements)]
        if group_two is None:
            g_one = methy_measurements[methy_measurements['name'].str.contains(group_one)]
            g_two = methy_measurements[methy_measurements['name'].str.contains(group_one) == False]
            pair = (g_one, g_two)
        else:
            g_one = methy_measurements[methy_measurements['annotation'].str.contains(group_one)]
            g_two = methy_measurements[methy_measurements['annotation'].str.contains(group_two)]
            pair = (g_one, g_two)
        return pair

    def to_list_of_dict(self, group):
        ret_val = []
        for ele in group:
            for i in self.methylation_types:

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
            partition_type = None
            group_one = "normal"
            group_two = "cancer"
            grouping = "all_pairs"
        return partition_type, group_one, group_two, grouping

    def compute(self, chromosome, start, end, additional=None):
        part_type, g_one, g_two, grouping = self.unpack_params(additional)

        methy_data = self.get_methy_data(start, end, chromosome)
        methy_corr_res = []

        group_one, group_two = self.partion(part_type, "")
        # exp_group_one = Methylation(start, end, chromosome, measurements=group_one.to_dict('records'))
        # group_one = [c for c in exp_group_one.columns if "_" in c]
        group_one = group_one['id'].tolist()
        group_one = self.to_list_of_dict(group_one)

        # Note in this case, part_type is always None, so this piece of code is probably not useful
        if part_type is not None:
            exp_group_two = Methylation(start, end, chromosome, measurements=group_two.to_dict('records'))
            group_two = [c for c in exp_group_two.columns if "_" in c]
            group_two = self.to_list_of_dict(group_two)
            group_pairs = [(x, y) for x in group_one for y in group_two]
        else:
            group_pairs = itertools.combinations(group_one, 2)
        # pvalue_list = []

        if not methy_data.empty:
            # loop through every possible combinations of methylation
            #for data_source_one, data_source_two in itertools.combinations(self.methylation_types, 2):
            for data_source_one, data_source_two in group_pairs:

                type1 = data_source_one["id"]
                type2 = data_source_two["id"]

                # check if there's data for these two methy types
                if type1 not in methy_data.columns or type2 not in methy_data.columns:
                    continue

                correlation_coefficient = pearsonr(methy_data[type1], methy_data[type2])
                data_range = {
                    'attr-one': [min(methy_data[type1]), max(methy_data[type1])],
                    'attr-two': [min(methy_data[type2]), max(methy_data[type2])]
                }
                attr = 'methylation' if self.methy_type == 'methy' else 'methylation diff'
                if correlation_coefficient[1] <= self.pval_threshold:
                    corr_obj = build_obj('correlation', attr,
                                        attr, True, data_source_one,
                                        data_source_two,
                                        correlation_coefficient[0],
                                        correlation_coefficient[1],
                                        ranges=data_range)
                    methy_corr_res.append(corr_obj)
            methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
                                    reverse=True)

            result = pd.Series(methy_corr_res)
            result = result.apply(pd.Series)

            parse_res = result

            # result = result.to_json(orient='records')
            # parse_res = json.loads(result)

        return parse_res
