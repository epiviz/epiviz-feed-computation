from .BaseStats import BaseStats

from itertools import combinations
from scipy.stats.stats import pearsonr

class Correlation(BaseStats):
    def __init__(self, measurements, pval_threshold):
        super(Correlation, self).__init__(measurements, pval_threshold)
        self.measurements = self.filter(measurements)

    def compute_stat(self, data1, data2, params=None):
        return pearsonr(data1, data2)

    def compute(self, chr, start, end, params):

        ## filter measurements first

        # group measurements - assuming all pairs here
        msets = combinations(self.measurements, 2)
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





