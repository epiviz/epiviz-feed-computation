import pandas as pd
import numpy as np
import math
import itertools
from scipy.stats.stats import pearsonr, ttest_ind
from old_feed.utils import build_obj
from old_feed.UI_functions import format_exp_methy_output
from stat_classes.StatMethod import StatMethod
from old_feed.data_functions import Methylation, Methylation_diff


class CorrelationMethy(StatMethod):

    def __init__(self, measurements, methy_type):
        super(). __init__(measurements)
        self.methylation_types = super().get_measurements_self(methy_type)
        self.methy_type = methy_type

    def get_methy_data(self, start, end, chromosome):
        if self.methy_type == "methy":
            methy_data = Methylation(start, end, chromosome, measurements=self.methylation_types)
        else:
            methy_data = Methylation_diff(start, end, chromosome, measurements=self.methylation_types)
        return methy_data

    def compute(self, chromosome, start, end):
        methy_data = self.get_methy_data(start, end, chromosome)
        methy_corr_res = []
        if not methy_data.empty:
            # loop through every possible combinations of methylation
            for data_source_one, data_source_two in itertools.combinations(self.methylation_types, 2):
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
                corr_obj = build_obj('correlation', 'methylation diff',
                                     'methylation diff', True, data_source_one,
                                     data_source_two,
                                     correlation_coefficient[0],
                                     correlation_coefficient[1],
                                     ranges=data_range)
                methy_corr_res.append(corr_obj)
            methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
                                    reverse=True)

            result = pd.Series(methy_corr_res)
            result = result.apply(pd.Series)

        return result
