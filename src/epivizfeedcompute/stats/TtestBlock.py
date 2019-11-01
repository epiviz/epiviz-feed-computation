import logging
import itertools
import pandas as pd
from scipy.stats import ttest_ind
from .BaseStats import BaseStats


class TtestBlock(BaseStats):

    def __init__(self, measurements, pval_threshold):
        super(TtestBlock, self).__init__(measurements, pval_threshold)
        self.measurements = measurements
        self.pval_threshold = pval_threshold

    def filter_measurements(self, params):
        '''
        Filters measurements for measurements with datatype needed by the current class

        Args:
            params (dict): holds a datatype field specifying the required datatype
        
        Returns:
            filtered measurments
        '''
        filtered = []

        for m in self.measurements:
            if m.annotation["datatype"] == "expression" or m.annotation["datatype"] == "peak":
                filtered.append(m)
        
        return filtered
    
    def get_transform_data(self, measurements, chr, start, end, params = None):
        '''
        Gets transform data
        Args: 
            measurements (list): list of measurements
        Returns:
            tuple of transformed data
        '''
        
        data = super(TtestBlock, self).get_transform_data(measurements, chr, start, end)
        block = []
        non_block = []
        out_block = data[0]
        for index, row in data[1].iterrows():
            in_block = data[0].where(((row.start <= data[0].start) & (data[0].start <= row.end)) | ((row.start <= data[0].end) & (data[0].end <= row.end))).dropna()

            out_curr_block = data[0].where(((row.start > data[0].start) | (data[0].start > row.end)) & ((row.start > data[0].end) | (data[0].end > row.end))).dropna()
            cols_to_use = curr_out.columns.difference(out_curr_block.columns)
            out_block = pd.merge(out_block, out_curr_block[cols_to_use], left_index=True, right_index=True)
            if len(in_block) > 0:
                block += in_block[measurements[0].mid].tolist()
        
        if len(out_block) > 0:
            non_block += out_block[measurements[0].mid].tolist()
        
            #       
            # for i, data_one_row in data[0].iterrows():    
            #     if (row.start <= data_one_row.start and data_one_row.start <= row.end) or (row.start <= data_one_row.end and data_one_row.end <= row.end):
            #         block.append(row[measurements[0].mid])
            #         in_block = True
            # if not in_block: 
            #     non_block.append(row[measurements[0].mid])
        return (block, non_block)

    def group_measurements(self, annotation):
        '''
        Groups measurement based on annotation if supplied. Otherwise preforms all pairs matching
        Args:
            annotation (str): specifies the grouping 
        Returns:
            all pairs between the two groups
        '''
        groups = {}
        
        if annotation == None:
            for m in self.measurements:
                if m.datatype in groups:
                    groups[m.datatype].append(m)
                else:
                    groups[m.datatype] = [m]
            #sorts items by datatype
            groups = sorted(groups.items(), key = lambda kv:kv[0])
            return itertools.product(groups[0][1], groups[1][1])
        else:
            for m in self.measurements:
                if m.annotation[annotation] in groups:
                    groups[annotation].append(m)
                else:
                    groups[annotation] = [m]
            return itertools.combinations(groups, len(groups.keys()))
    
    def compute_stat(self, data1, data2, params=None):
        '''
        Stat method

        Args:
            data1 (dict): data in group 1
            data2 (dict): data in group 2
        Returns:
            ttest value and pvalue
        '''
        
        return ttest_ind(data1, data2, equal_var=False)

    def compute(self, chr, start, end, params=None):
        '''
        
        Computes statistical method on the given measurement
        Args:
            chr (str): chromosome
            start (int): genomic start
            end (int): genomic end
            params (dict): additional params contains annotations and datatype field
        Returns:
            dataframe of results
            
        '''
        self.measurements = self.filter_measurements(params)
        msets = self.group_measurements(params["annotation"])
        results = []
        
        for (m1, m2) in msets:
            data1, data2 = self.get_transform_data([m1, m2], chr, start, end)
            value, pvalue = self.compute_stat(data1, data2)
            if pvalue <= self.pval_threshold:
                results.append(
                    {
                        'measurements': (m1, m2),
                        'test': 'ttest',
                        'value': value,
                        'pvalue': pvalue
                    }
                )

        sorted_results = sorted(results, key=lambda x: x['value'], reverse=True)

        return self.toDataFrame(sorted_results)
