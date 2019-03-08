import pandas as pd
import numpy as np
import math
import itertools
from scipy.stats.stats import pearsonr, ttest_ind
from old_feed.utils import build_obj
from old_feed.UI_functions import format_exp_methy_output
from stat_classes.stat_method import StatMethod
from old_feed.data_functions import Gene_data, Methylation


class CorrelationExp(StatMethod):

    def __init__(self, measurements):
        super(). __init__(measurements)
        self.gene_types = super().get_measurements_self("gene")

    def compute(self, chromosome, start, end):
        expression_data = Gene_data(start, end, chromosome, measurements=self.gene_types)
        corr_list = []
        # pvalue_list = []
        for data_source_one, data_source_two in itertools.combinations(self.gene_types, 2):
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

        return corr_res
