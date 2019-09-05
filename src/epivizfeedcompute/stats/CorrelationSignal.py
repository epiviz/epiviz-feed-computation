from .BaseStats import BaseStats

from itertools import combinations
from scipy.stats.stats import pearsonr

class CorrelationSignal(BaseStats):
    def __init__(self, measurements, pval_threshold):
        super(CorrelationSignal, self). __init__(measurements, pval_threshold)
        self.measurements = self.filter(measurements)

    def compute(self, chr, start, end, params):

        ## filter measurements first

        # group measurements - assuming all pairs here
        msets = combinations(self.measurements, 2)
        results = []

        for (m1, m2) in msets:
            corr, pvalue = pearsonr(methy_data[type1], methy_data[type2])

            if pval <= pval_threshold:
                results.append(
                    {
                        'measurements': (m1, m2),
                        'test': 'correlation',
                        'value': corr,
                        'pvalue': pval
                    }
                )

        sorted_results = sorted(results, key=lambda x: x['value'], reverse=True)

        return self.toDataFrame(sorted_results)





