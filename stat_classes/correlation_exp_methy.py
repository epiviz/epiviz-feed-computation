import pandas as pd
import numpy as np
import math
import itertools
from scipy.stats.stats import pearsonr
from utils import build_exp_methy_obj
from UI_functions import format_exp_methy_output
from stat_classes.stat_method import stat_method
from data_functions import Gene_data, Methylation, Methylation_diff


class CorrelationExpMethy(stat_method):

    def __init__(self, measurements, methy_name):
        super(). __init__(measurements)
        self.datasource_gene_types = super().get_measurements_self("gene")
        self.datasource_methy_types = super().get_measurements_self(methy_name)
        self.methy_name = methy_name

    def find_expression_block(self, regions_of_interest, methy_data):
        mean = pd.DataFrame(columns=[methy_type["id"] for methy_type in self.datasource_methy_types])
        for index, row in regions_of_interest.iterrows():
            mean = mean.append(methy_data.where(((row.start <= methy_data.start) & (methy_data.start <= row.end)) | ((row.start <= methy_data.end) & (methy_data.end <= row.end))).mean().fillna(0), ignore_index=True)
        return mean

    def correlation_calc(self, exp_data, methy_mean, tissue_pair):
        ret_obj = []
        # use the difference of types (normal and tumor) for the same tissue
        # this need to be changed to accomadate multiple attributes
        exp_type1 = tissue_pair[0]
        exp_type2 = tissue_pair[1]

        tissue_type = exp_type1.split("___")[0]
        # for datasource_methy_type in datasource_methy_types:
        methy_type_norm = tissue_type + "_" + exp_type1.split("___")[1]
        methy_type_canc = tissue_type + "_cancer"
        expression_diff = np.subtract(exp_data[exp_type1], exp_data[exp_type2])
        if self.methy_name == "methy":
            correlation_coefficient_norm = pearsonr(methy_mean[methy_type_norm], expression_diff)
            correlation_coefficient_canc = pearsonr(methy_mean[methy_type_canc], expression_diff)
            if not math.isnan(correlation_coefficient_norm[0]):
                # format the data into list of json objects for plots
                data = format_exp_methy_output(expression_diff, methy_mean[methy_type_norm], "Expression diff " + tissue_type,  'Collapsed Methylation Diff ' + tissue_type)

                data_range = {
                    'attr-one': [min(expression_diff), max(expression_diff)],
                    'attr-two': [min(methy_mean[methy_type_norm]), max(methy_mean[methy_type_norm])]
                }

                corr_obj = build_exp_methy_obj('correlation', 'expression diff',
                                               'methylation diff', True, "Expression diff " + tissue_type,
                                               'Collapsed Methylation Diff ' + tissue_type,
                                               correlation_coefficient_norm[0],
                                               correlation_coefficient_norm[1],
                                               data=data, ranges=data_range)
                data_canc = format_exp_methy_output(expression_diff, methy_mean[methy_type_canc], "Expression diff " + tissue_type,  'Collapsed Methylation Diff ' + tissue_type)

                data_canc_range = {
                    'attr-one': [min(expression_diff), max(expression_diff)],
                    'attr-two': [min(methy_mean[methy_type_canc]), max(methy_mean[methy_type_canc])]
                }

                corr_obj_cancer = build_exp_methy_obj('correlation', 'expression diff',
                                               'methylation diff', True, "Expression diff " + tissue_type + "cancer",
                                               'Collapsed Methylation Diff ' + tissue_type,
                                               correlation_coefficient_canc[0],
                                               correlation_coefficient_canc[1],
                                               data=data_canc, ranges=data_range)
                ret_obj = [corr_obj, corr_obj_cancer]
        else:
                correlation_coefficient = pearsonr(methy_mean[tissue_type], expression_diff)
                data_out = format_exp_methy_output(expression_diff, methy_mean[tissue_type], "Expression diff " + tissue_type,  'Collapsed Methylation Diff ' + tissue_type)

                data_range = {
                    'attr-one': [min(expression_diff), max(expression_diff)],
                    'attr-two': [min(methy_mean[tissue_type]), max(methy_mean[tissue_type])]
                }

                corr_obj = build_exp_methy_obj('correlation', 'expression diff',
                                               'methylation diff', True, "Expression diff " + tissue_type,
                                               'Collapsed Methylation Diff ' + tissue_type,
                                               correlation_coefficient[0],
                                               correlation_coefficient[1],
                                               data=data_out, ranges=data_range)
                ret_obj = [corr_obj]

        return ret_obj

    def get_methy_data(self, chromosome, start, end):
        if self.methy_name == "methy":
            methy_data = Methylation(start, end, chromosome, measurements=self.datasource_methy_types)

        else:
            methy_data = Methylation_diff(start, end, chromosome, measurements=self.datasource_methy_types)
        return methy_data

    def compute(self, start, end, chromosome, downstream=3000, upstream=1000):
        exp_data = Gene_data(start, end, chromosome, measurements=self.datasource_gene_types)
        methy_data = self.get_methy_data(chromosome, start, end)
        results = []
        exp_regions = np.array([exp_data.start, exp_data.end]).transpose()
        down_up = np.ones(exp_regions.shape) * np.array([downstream, -1 * upstream])

        regions_of_interest = pd.DataFrame(np.subtract(exp_regions, down_up), columns=["start", "end"])
        print(regions_of_interest)
        #names cols to methy types in datasource methy types
        methy_mean = self.find_expression_block(regions_of_interest, methy_data)

        # sorts gene_types by tissue type
        sorted_gene_types = sorted(self.datasource_gene_types, key=lambda x: x["id"])
        # gets list of dicts containing the same info as gene type
        tissue_types = [item["id"] for item in sorted_gene_types]
        tissue_types = [tissue_types[x:x+2] for x in range(0, len(tissue_types), 2)]
        for tissue_pair in tissue_types:
            corr_obj_res = self.correlation_calc(exp_data, methy_mean, tissue_pair)
            if self.methy_name == "methy":
                results.append(corr_obj_res[0])
                results.append(corr_obj_res[1])
            else:
                results.append(corr_obj_res[0])

        # corr_res = sorted(corr_res, key=lambda x: x['value'], reverse=True)
        # corr_cancer = sorted(corr_res_cancer, key=lambda x: x['value'], reverse=True)
        corr_result = pd.Series(results)
        corr_result = corr_result.apply(pd.Series)

        print(corr_result)
        return corr_result


