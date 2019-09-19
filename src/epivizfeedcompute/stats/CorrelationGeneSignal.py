from .BaseStats import BaseStats

from itertools import combinations
from scipy.stats.stats import pearsonr

class CorrelationGeneSignal(BaseStats):
    def __init__(self, measurements, pval_threshold):
        super(Correlation, self).__init__(measurements, pval_threshold)
        self.measurements = measurements
    def filter_measurements(self, params):
        filtered = []
        for m in self.measurements:
            if m.datatype == "expression" or m.datatype == "signal":
                filtered.append(m)
        return filtered
    
    def get_transform_data(self, measurements):
        data = super(Correlation, self).get_transform_data(measurements)

        mean = []
        for index, row in data[0].iterrows():
            mean = mean.append(data[1].where(((row.start - upstream <= data[1].start) & (data[1].start <= row.end + downstream)) | ((row.start - upstream <= data[1].end) & (data[1].end <= row.end + downstream))).mean().fillna(0), ignore_index=True)
        return (data[0], mean)
        # data[1] = data[1].set_index(['start', 'end'])
        # data[1].index = pd.IntervalIndex.from_tuples(data[1].index)

        # # assuming data[0] is expression and data[1] is methylation
        # for index, row in data[0].iterrows():

        
    def compute(self, chr, start, end, params, upstream=1000, downstream=3000):
        self.measurements = self.filter(self.measurements)
        self.upstream = upstream
        self.downstream = downstream
        return super(Correlation, self).compute(chr, start, end, params)





