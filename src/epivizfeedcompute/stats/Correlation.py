from .BaseStats import BaseStats

from itertools import combinations
from scipy.stats.stats import pearsonr

class Correlation(BaseStats):
    def __init__(self, measurements, pval_threshold):
        super(Correlation, self).__init__(measurements, pval_threshold)
        self.measurements = measurements

    def filter_measurements(self, params):
        filtered = []
        for m in self.measurements:
            if m.datatype == params.datatype:
                filtered.append(m)
        return filtered
        
    def compute_stat(self, data1, data2, params=None):
        return pearsonr(data1, data2)

    def group_measurements(self, annotation):
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
            
    def compute(self, chr, start, end, params):

        ## filter measurements first
        self.measurements = self.filter_measurements(params)
        # group measurements - assuming all pairs here
        msets = self.group_measurements(params.annotation)
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





