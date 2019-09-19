import numpy as np
import pandas as pd
import logging
import math
import json
import itertools
from epivizFeed.utils import build_obj
from scipy.stats import fisher_exact
from epivizFeed.StatMethod import StatMethod
from epivizFeed.data_functions import Block_data
from .BaseStats import BaseStats


class OverlapBlock(BaseStats):

    def __init__(self, measurements):
        super(OverlapBlock, self).__init__(measurements)
        self.measurements = measurements
    def filter_measurements(self, params):
        filtered = []
        for m in self.measurements:
            if m.datatype == params.datatype:
                filtered.append(m)
        return filtered
        
    def get_transform_data(self, measurements, params=None):
        data = super(OverlapBlock, self).get_transform_data(measurements)
        start = params["start"]
        end = params["end"]
        b_one_ind = 0
        b_two_ind = 0
        d1 = (0,0)
        d2 = (0,0)
        d1_sum = 0
        d2_sum = 0
        while b_one_ind < len(data[0]) and b_two_ind < len(data[1]):

            tissue_one_start = max(float(start), data[0]['start'][b_one_ind])
            tissue_two_start = max(float(start), data[1]['start'][b_two_ind])
            tissue_one_end = min(float(end), data[0]['end'][b_one_ind])
            tissue_two_end = min(float(end), data[1]['end'][b_two_ind])
            # there is an overlap
            if tissue_one_end < tissue_two_start or tissue_two_end < tissue_one_start:
                d1[1] += tissue_one_end - tissue_one_start
                d2[1] += tissue_two_end - tissue_two_start
            else:
                common_end = min(tissue_two_end, tissue_one_end)
                common_start = max(tissue_one_start, tissue_two_start)
                d1[0] += common_end - common_start
                d2[0] += common_end - common_start
                if tissue_one_start < tissue_two_start < tissue_two_end < tissue_one_end:
                    d1[1] += tissue_two_start - tissue_one_start  + tissue_one_end - tissue_two_end 
                elif tissue_two_start < tissue_one_start < tissue_one_end < tissue_two_end:
                    d2[1] += tissue_one_start - tissue_two_start  + tissue_two_end - tissue_one_end                        
                elif tissue_one_end > tissue_two_start:
                    d1[1] += tissue_two_start - tissue_one_start 
                elif tissue_two_end > tissue_one_start
                    d2[1] += tissue_one_start - tissue_two_start 
            # block tissue two is larger
            if tissue_two_start >= tissue_one_end or tissue_two_end > tissue_one_end:
                b_one_ind += 1
                d1_sum += tissue_one_end - tissue_one_start
            
            # block tissue one is larger
            if tissue_one_start >= tissue_two_end or tissue_two_end <= tissue_one_end:
                b_two_ind += 1
                d2_sum += tissue_two_end - tissue_two_start

        return ([d1[0], d1[1]/d1_sum], [d2[0], d2[1]/d2_sum])
   
    def compute_stat(self, data1, data2, params=None):
        return fisher_exact(np.array([data1, data2]))

    def compute(self, chr, start, end, params):
        self.measurements = self.filter(params)
        msets = combinations(self.measurements, 2)
        results = []
        
        for (m1, m2) in msets:
            data1, data2 = self.get_transform_data([m1, m2])
            value, pvalue = self.compute_stat(data1, data2)
            if pvalue <= self.pval_threshold:
                results.append(
                    {
                        'measurements': (m1, m2),
                        'test': 'overlap',
                        'value': value,
                        'pvalue': pvalue
                    }
                )

        sorted_results = sorted(results, key=lambda x: x['value'], reverse=True)

        return self.toDataFrame(sorted_results)
