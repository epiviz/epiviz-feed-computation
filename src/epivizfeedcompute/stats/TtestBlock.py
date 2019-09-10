import logging
from scipy.stats import ttest_ind
from .BaseStats import BaseStats


class TtestBlock(BaseStats):

    def __init__(self, measurements):
        super(TtestBlock, self).__init__(measurements)
        self.measurements = self.filter(measurements)

     def get_transform_data(self, measurements):
        data = super(Correlation, self).get_transform_data(measurements)
        block = []
        non_block = []
        for index, row in data[0].iterrows():
            if data[1].where(((row.start <= data[1].start) & (data[1].start <= row.end)) | ((row.start <= data[1].end) & (data[1].end <= row.end)):
                block.append(row[measurements[0].mid])
            else:
                non_block.append(row[measurements[0].mid])
        return (block, non_block)
   
    def compute_stat(self, data1, data2, params=None):
        
        return ttest_ind(data1, data2, equal_var=False)

    def compute(self, chr, start, end, params):
        msets = combinations(self.measurements, 2)
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
