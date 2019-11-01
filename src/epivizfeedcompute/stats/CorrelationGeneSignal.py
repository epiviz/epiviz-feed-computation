import pandas as pd
from .BaseStats import BaseStats
import itertools
from scipy.stats.stats import pearsonr
from .Correlation import Correlation

class CorrelationGeneSignal(Correlation):
    def __init__(self, measurements, pval_threshold):
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
            if m.annotation["datatype"] == "expression" or m.annotation["datatype"] == "signal":
                filtered.append(m)
        return filtered
    
    def get_transform_data(self, measurements, chr, start, end, params = None):
        '''
        Gets transform data
        Args: 
            measurements (list): list of measurements
        Returns:
            tuple transformed data
        '''
        data = []
        mids = []
        mean = []

        for m in measurements:
            data.append(m.get_data(chr, start, end)[0])
            mids.append(m.mid)
        for index, row in data[0].iterrows():
            mean.append(data[1].where(((row.start - self.upstream <= data[1].start) & (data[1].start <= row.end + self.downstream)) | ((row.start - self.upstream <= data[1].end) & (data[1].end <= row.end + self.downstream))).mean().fillna(0))
        return (data[0][mids[0]], pd.DataFrame(mean)[mids[1]])
        
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
    
    def compute(self, chr, start, end, params= None, upstream=1000, downstream=3000):
        '''
        Computes statistical method on the given measurement
        
        Args:
            chr (str): chromosome
            start (int): genomic start
            end (int): genomic end
            params (dict): additional params contains annotations and datatype field
            upstream (int): number of base pairs upstream
            downstream (int): number of base pairs downstream
        Returns:
            dataframe of results
            
        '''
        self.measurements = self.filter_measurements(self.measurements)
        self.upstream = upstream
        self.downstream = downstream
        return super().compute(chr, start, end, params)





