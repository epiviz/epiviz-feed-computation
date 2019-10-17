from .BaseStats import BaseStats

import itertools 
from scipy.stats.stats import pearsonr

class Correlation(BaseStats):
    def __init__(self, measurements, pval_threshold):
        '''
        Correlation class initialization
        Args: 
            measurements (list): list of measurements
            pval_threshold (double): p value threshold
        '''
        super(Correlation, self).__init__(measurements, pval_threshold)
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
            if m["datatype"] == params["datatype"]:
                filtered.append(m)
        return filtered
        
    def compute_stat(self, data1, data2, params=None):
        '''
        Stat method

        Args:
            data1 (dict): data in group 1
            data2 (dict): data in group 2
        Returns:
            pearsonr coefficent and pvalue
        '''
        
        return pearsonr(data1, data2)

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
            return itertools.combinations(self.measurements, 2)
        else:
            for m in self.measurements:
                if m.annotation[annotation] in groups:
                    groups[annotation].append(m)
                else:
                    groups[annotation] = [m]
            return itertools.combinations(groups, len(groups.keys()))
            
    def compute(self, chr, start, end, params):
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
            data1, data2 = self.get_transform_data([m1, m2])
            corr, pvalue = self.compute_stat(data1, data2)

            if pvalue <= self.pval_threshold:
                results.append(
                    {
                        'measurements': (m1, m2),
                        'test': 'correlation',
                        'value': corr,
                        'pvalue': pvalue
                    }
                )

        sorted_results = sorted(results, key=lambda x: x['value'], reverse=True)

        return self.toDataFrame(sorted_results)





