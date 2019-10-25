import pandas as pd

class BaseStats(object):
    '''
    Base class for all stat methods
    '''
    def __init__(self, measurements, pval_threshold):
        '''
        Base class initialization
        Args: 
            measurements (list): list of measurements
            pval_threshold (double): p value threshold
        '''
        self.measurements = measurements
        self.pval_threshold = pval_threshold

    def filter_measurements(self):
        raise Exception("NotImplementedException")

    def group_measurements(self, annotation):
        raise Exception("NotImplementedException")

    def compute(self, chr, start, end, params=None):
        ## steps here
        # 1. GroupBy
        # 2. getData for measurements in each group
        # 3. Calculate Test
        raise Exception("NotImplementedException")

    def get_transform_data(self,measurements, chr, start, end, params = None):
        '''
        Gets transform data
        Args: 
            measurements (list): list of measurements
        Returns:
            tuple transformed data
        '''
        data = []

        for m in measurements:
            print(m)
            data.append(m.get_data(chr, start, end)[0])
        return tuple(data)

    def toDataFrame(self, results):
        '''
        Converts results to dataframe
        Args: 
            measurements (list): list of measurements
            pval_threshold (double): p value threshold
        Returns:
            Dataframe of result data
        '''
        measurement1 = [r['measurements'][0].name for r in results] 
        measurement2 = [r['measurements'][1].name for r in results] 
        value = [r['value'] for r in results] 
        pvalue = [r['pvalue'] for r in results] 
        return pd.DataFrame(
            data = {
                "measurement1" : measurement1,
                "measurement2" : measurement2,
                "value" : value,
                "pvalue" : pvalue
            }
        )



    