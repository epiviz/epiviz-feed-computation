import logging
import math
import itertools
from scipy.stats import fisher_exact
from .BaseStats import BaseStats

class OverlapBlock(BaseStats):

    def __init__(self, measurements):
        super(OverlapBlock, self).__init__(measurements)
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
            if m.datatype == params.datatype:
                filtered.append(m)
        return filtered
        
    def get_transform_data(self, measurements, params=None):
        '''
        Gets transform data
        Args: 
            measurements (list): list of measurements
        Returns:
            tuple of transformed data
        '''
        
        data = super(OverlapBlock, self).get_transform_data(measurements)
        start = params["start"]
        end = params["end"]
        b_one_ind = 0
        b_two_ind = 0
        data_one = (0,0)
        data_two = (0,0)
        data_one_sum = 0
        data_two_sum = 0
        while b_one_ind < len(data[0]) and b_two_ind < len(data[1]):
            block_one = (max(float(start), data[0]['start'][b_one_ind]), min(float(end), data[0]['end'][b_one_ind]))
            block_two = (max(float(start), data[1]['start'][b_two_ind]), min(float(end), data[1]['end'][b_two_ind]))
            # there is no overlap
            if block_one[1] < block_two[0] or block_two[1] < block_one[0]:
                data_one[1] += block_one[1] - block_one[0]
                data_two[1] += block_two[1] - block_two[0]
            else:
                data_one[0] += min(block_two[1], block_one[1]) - max(block_one[0], block_two[0])
                data_two[0] += min(block_two[1], block_one[1]) - max(block_one[0], block_two[0])

                if block_one[0] < block_two[0] < block_two[1] < block_one[1]:
                    data_one[1] += block_two[0] - block_one[0]  + block_one[1] - block_two[1] 

                elif block_two[0] < block_one[0] < block_one[1] < block_two[1]:
                    data_two[1] += block_one[0] - block_two[0]  + block_two[1] - block_one[1]     

                elif block_one[1] > block_two[0]:
                    data_one[1] += block_two[0] - block_one[0] 
                else 
                    data_two[1] += block_one[0] - block_two[0] 
            # block tissue two is larger
            if block_two[0] >= block_one[1] or block_two[1] > block_one[1]:
                b_one_ind += 1
                data_one_sum += block_one[1] - block_one[0]
            else:
                b_two_ind += 1
                data_two_sum += block_two[1] - block_two[0]

        return ([data_one[0], data_one[1]/data_one_sum], [data_two[0], data_two[1]/data_two_sum])
   
    def compute_stat(self, data1, data2, params=None):
        '''
        Stat method

        Args:
        redo
            data1 (dict): holds a datatype field specifying the required datatype
            data2 (dict): holds a datatype field specifying the required datatype
        Returns:
            oddsratio and pvalue
        '''
        
        return fisher_exact(np.array([data1, data2]))

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
                        'test': 'overlap',
                        'value': value,
                        'pvalue': pvalue
                    }
                )

        sorted_results = sorted(results, key=lambda x: x['value'], reverse=True)

        return self.toDataFrame(sorted_results)
