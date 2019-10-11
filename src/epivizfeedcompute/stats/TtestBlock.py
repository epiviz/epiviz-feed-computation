import logging
from scipy.stats import ttest_ind
from .BaseStats import BaseStats


class TtestBlock(BaseStats):

    def __init__(self, measurements):
        super(TtestBlock, self).__init__(measurements)
        self.measurements = measurements
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
            if m.datatype == "expression" or m.datatype == "peak":
                filtered.append(m)
        return filtered
    
     def get_transform_data(self, measurements):
        '''
        Gets transform data
        Args: 
            measurements (list): list of measurements
        Returns:
            tuple of transformed data
        '''
        
        data = super(TtestBlock, self).get_transform_data(measurements)
        block = []
        non_block = []
        for index, row in data[0].iterrows():
            if data[1].where(((row.start <= data[1].start) & (data[1].start <= row.end)) | ((row.start <= data[1].end) & (data[1].end <= row.end)):
                block.append(row[measurements[0].mid])
            else:
                non_block.append(row[measurements[0].mid])
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
            return combinations(self.measurements, 2)
        else:
            for m in self.measurements:
                if m.annotation[annotation] in groups:
                    groups[annotation].append(m)
                else:
                    groups[annotation] = [m]
            return combinations(groups, len(groups.keys()))
    
    def compute_stat(self, data1, data2, params=None):
        '''
        Stat method

        Args:
        redo
            data1 (dict): holds a datatype field specifying the required datatype
            data2 (dict): holds a datatype field specifying the required datatype
        Returns:
            ttest value and pvalue
        '''
        
        return ttest_ind(data1, data2, equal_var=False)

    def compute(self, chr, start, end, params):
        '''
        Computes stat method

        Args:
            chr (str): chromosome
            start (int): genomic start
            end (int): genomic end
            params (dict): additional params contains annotations and datatype field
        Returns:
            dataframe of results
            
        '''
        self.measurements = self.filter(params)
        msets = self.group_measurements(params.annotation)
        results = []
        
        for (m1, m2) in msets:
            data1, data2 = self.get_transform_data([m1, m2])
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
