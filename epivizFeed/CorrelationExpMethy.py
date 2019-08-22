import pandas as pd
import numpy as np
import math
import json
import itertools
from scipy.stats.stats import pearsonr
from old_feed.utils import build_exp_methy_obj
from old_feed.UI_functions import format_exp_methy_output
from epivizFeed.StatMethod import StatMethod
from old_feed.data_functions import Gene_data, Methylation, Methylation_diff


class CorrelationExpMethy(StatMethod):

    def __init__(self, measurements, methy_name, pval_threshold):
        super(CorrelationExpMethy, self). __init__(measurements, pval_threshold)
        self.datasource_gene_types = super(CorrelationExpMethy, self).get_measurements_self("gene")
        self.datasource_methy_types = super(CorrelationExpMethy, self).get_measurements_self(methy_name)
        self.methy_name = methy_name

    def to_list_of_dict(self, ele):
        ret_val = None
        for i in self.datasource_gene_types:
            if i["id"] == ele:
                ret_val = i
        return ret_val

    def partion(self, type, group_one, group_two=None):
        #attributes other than tissues must be specifed in form [attr_name, (attr_val 1, attr_val 2)]
        # where attributes values must be distinct strings
        pair = None
        m = pd.DataFrame()
        for measurement in self.measurements:
            m = m.append(measurement, ignore_index=True)
        exp_measurements = m[m['name'].str.contains('Expression')]

        if group_two is None:
            g_one = exp_measurements[exp_measurements['name'].str.contains(group_one)]
            g_two = exp_measurements[exp_measurements['name'].str.contains(group_one) == False]
            pair = (g_one, g_two)
        else:
            g_one = exp_measurements[exp_measurements['annotation'].str.contains(group_one)]
            g_two = exp_measurements[exp_measurements['annotation'].str.contains(group_two)]
            pair = (g_one, g_two)
        return pair

    def grouping(self, group_one, group_two, all_pairs=True):
        group_pairs = []

        if all_pairs:
            g_one = [c for c in group_one.columns if "_" in c]
            g_two = [c for c in group_two.columns if "_" in c]

            group_pairs = [(x, y) for x in g_one for y in g_two]
            group_pairs = [(self.to_list_of_dict(x), self.to_list_of_dict(y)) for x, y in group_pairs]
        else:
            a = group_one.columns
            b = group_two.columns

            group_pairs = [(x, y) for x, y in zip(a, b) if x not in b and y not in a]
            group_pairs = [(self.to_list_of_dict(x), self.to_list_of_dict(y)) for x, y in group_pairs]

        return group_pairs

    def find_expression_block(self, regions_of_interest, methy_data):
        mean = pd.DataFrame(columns=[methy_type["id"] for methy_type in self.datasource_methy_types])
        for index, row in regions_of_interest.iterrows():
            mean = mean.append(methy_data.where(((row.start <= methy_data.start) & (methy_data.start <= row.end)) | ((row.start <= methy_data.end) & (methy_data.end <= row.end))).mean().fillna(0), ignore_index=True)
        return mean

    def create_corr_diff_obj(self, exp_diff, methy_col, tissue_type, corr_coef):
        data = format_exp_methy_output(exp_diff, methy_col, "Expression diff "
            + tissue_type,  'Collapsed Methylation Diff ' + tissue_type)

        data_range = {
            'attr-one': [min(exp_diff), max(exp_diff)],
            'attr-two': [min(methy_col), max(methy_col)]
        }
        corr_obj = build_exp_methy_obj('correlation', 'expression diff',
                                       'methylation diff', True, "Expression diff " + tissue_type,
                                       'Collapsed Methylation Diff ' + tissue_type,
                                       corr_coef[0], corr_coef[1],
                                       data=data, ranges=data_range)
        return corr_obj

    def create_corr_obj(self, expression, methy_mean, exp_type, methy_type, corr_coef):
        # format the data into list of json objects for plots
        data = format_exp_methy_output(expression, methy_mean, "Expression " + exp_type.split('___')[0] + '_' + exp_type.split('___')[1], 
        "Average Probe level Meth " + methy_type)

        data_range = {
            'attr-one': [min(expression),
                        max(expression)],
            'attr-two': [min(methy_mean),
                        max(methy_mean)]
        }
        # if correlation_coefficient[1] <= self.pval_threshold:
        corr_obj = build_exp_methy_obj('correlation', 'expression',
                                    'methylation', True, "Expression " + exp_type.split('___')[0] + '_' + exp_type.split('___')[1], 
                                    "Average Probe level Meth " + methy_type,
                                    corr_coef[0],
                                    corr_coef[1],
                                    data=data, ranges=data_range)
        return corr_obj

    def unpack_params(self, additional):
        if additional is not None:
            partition_type = additional["partition_type"]
            group_one = additional["group_one"]
            group_two = additional["group_two"]
            grouping = additional["grouping"]
        else:
            partition_type = "condition"
            group_one = "normal"
            group_two = "cancer"
            grouping = "all_pairs"
        return partition_type, group_one, group_two, grouping

    def correlation_calc(self, exp_data, methy_mean, tissue_pair):
        ret_obj = []
        # use the difference of types (normal and tumor) for the same tissue
        # this need to be changed to accomadate multiple attributes
        exp_type1 = tissue_pair[0]["id"]
        exp_type2 = tissue_pair[1]["id"]
        tissue_type = exp_type1.split("___")[0]
        # for datasource_methy_type in datasource_methy_types:
        methy_type_norm = tissue_type + "_" + exp_type1.split("___")[1]
        methy_type_canc = tissue_type + "_cancer"
        expression_diff = np.subtract(exp_data[exp_type1], exp_data[exp_type2])
        res_obj = []

        if self.methy_name == "methy":
            correlation_coefficient_norm = pearsonr(methy_mean[methy_type_norm], exp_data[exp_type1])
            correlation_coefficient_canc = pearsonr(methy_mean[methy_type_canc], exp_data[exp_type2])
            correlation_coefficient_3 = pearsonr(methy_mean[methy_type_norm], exp_data[exp_type2])
            correlation_coefficient_4 = pearsonr(methy_mean[methy_type_canc], exp_data[exp_type1])
            if not math.isnan(correlation_coefficient_norm[0]) and correlation_coefficient_norm[1] <= self.pval_threshold:
                # format the data into list of json objects for plots
                corr_obj_1 = self.create_corr_obj(exp_data, methy_mean[methy_type_norm], exp_type1, methy_type_norm, correlation_coefficient_norm)
                ret_obj.append(corr_obj_1)
            if not math.isnan(correlation_coefficient_canc[1]) and correlation_coefficient_canc[1] <= self.pval_threshold:
                corr_obj_2 = self.create_corr_obj(exp_data, methy_mean[methy_type_canc], exp_type2, methy_type_canc, correlation_coefficient_canc)
                ret_obj.append(corr_obj_2)
            if not math.isnan(correlation_coefficient_3[1]) and correlation_coefficient_3[1] <= self.pval_threshold:
                corr_obj_3 = self.create_corr_obj(exp_data, methy_mean[methy_type_norm], exp_type2, methy_type_norm, correlation_coefficient_3)
                ret_obj.append(corr_obj_3)
            if not math.isnan(correlation_coefficient_4[1]) and correlation_coefficient_4[1] <= self.pval_threshold:
                corr_obj_4 = self.create_corr_obj(exp_data, methy_mean[methy_type_canc], exp_type1, methy_type_canc, correlation_coefficient_4)
                ret_obj.append(corr_obj_4)
        else:
            correlation_coefficient = pearsonr(methy_mean[tissue_type], expression_diff)
            if not math.isnan(correlation_coefficient[0]) and correlation_coefficient[1] <= self.pval_threshold:
                # format the data into list of json objects for plots
                corr_obj = self.create_corr_diff_obj(expression_diff, methy_mean[tissue_type], tissue_type, correlation_coefficient)
                ret_obj = [corr_obj]

        return ret_obj

    def get_methy_data(self, chromosome, start, end):
        if self.methy_name == "methy":
            methy_data = Methylation(start, end, chromosome, measurements=self.datasource_methy_types)

        else:
            methy_data = Methylation_diff(start, end, chromosome, measurements=self.datasource_methy_types)
        return methy_data

    def same_tissue_match(self, pair):
        if pair[0]['id'].split('___')[0] == pair[1]['id'].split('___')[0]:
            return True
        else:
            return False    

    def compute(self, chromosome, start, end, additional=None, downstream=3000, upstream=1000):
        part_type, g_one, g_two, grouping = self.unpack_params(additional)
        if grouping == "all_pairs":
            grouping = True
        else:
            grouping = False

        exp_data = Gene_data(start, end, chromosome, measurements=self.datasource_gene_types)
        methy_data = self.get_methy_data(chromosome, start, end)
        group_one, group_two = self.partion(part_type, g_one)
        exp_group_one = Gene_data(start, end, chromosome, measurements=group_one.to_dict('records'))
        exp_group_two = Gene_data(start, end, chromosome, measurements=group_two.to_dict('records'))

        group_pairs = self.grouping(exp_group_one, exp_group_two, grouping)
        # if self.methy_name == 'methy_diff':
        #filter group_pairs with only the same tissue
        group_pairs = filter(self.same_tissue_match, group_pairs)
        results = []
        exp_regions = np.array([exp_data.start, exp_data.end]).transpose()
        down_up = np.ones(exp_regions.shape) * np.array([downstream, -1 * upstream])

        regions_of_interest = pd.DataFrame(np.subtract(exp_regions, down_up), columns=["start", "end"])
        #names cols to methy types in datasource methy types
        methy_mean = self.find_expression_block(regions_of_interest, methy_data)

        # # sorts gene_types by tissue type
        # sorted_gene_types = sorted(self.datasource_gene_types, key=lambda x: x["id"])
        # # gets list of dicts containing the same info as gene type
        # tissue_types = [item["id"] for item in sorted_gene_types]
        # tissue_types = [tissue_types[x:x+2] for x in range(0, len(tissue_types), 2)]

        for tissue_pair in group_pairs:
            corr_obj_res = self.correlation_calc(exp_data, methy_mean, tissue_pair)
            # if self.methy_name == "methy":
            results.extend(corr_obj_res)
            # else:
            #     results.extend(corr_obj_res)

        results = sorted(results, key=lambda x: abs(x['value']),
                                reverse=True)
        corr_result = pd.Series(results)
        corr_result = corr_result.apply(pd.Series)

        parse_res = corr_result
        # corr_result = corr_result.to_json(orient='records')
        # parse_res = json.loads(corr_result)

        return parse_res


