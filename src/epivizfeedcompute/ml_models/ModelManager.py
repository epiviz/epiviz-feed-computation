import tensorflow as tf
import keras
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Activation, Flatten, Conv1D, MaxPooling1D
from keras.layers.core import Reshape
from .BaseModel import BaseModel

class ModelManager():

    def __init__ (self, measurements, pval_threshold, MeasurementTypes=None):
        self.models = {}
        self.preprocess_functions = {}
        self.pval_threshold = pval_threshold
        self.measurements = measurements
        self.expectedMeasurement = self.query_measurements(measurements, "lung___tumor")

    def query_measurements(self, measurements, query=None):
        if query == None:
            return measurements
        results = []
        for m in measurements:
            if m.mid == query:
                results.append(m)
        return results

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
        measurement1Data = [r['measurements'][0] for r in results] 
        measurement2Data = [r['measurements'][1] for r in results] 
        value = [r['value'] for r in results] 
        pvalue = [r['pvalue'] for r in results] 
        test = [r['test'] for r in results]
        return pd.DataFrame(
            data = {
                "test":test,
                "measurement1" : measurement1,
                "measurement2" : measurement2,
                
                "measurement1Data" : measurement1Data,
                "measurement2Data" : measurement2Data,
                "value" : value,
                "pvalue" : pvalue
            }
        )

    def Add(self, modelName, model, preprocess_function):
        self.models[modelName] = BaseModel(self.measurements, model, preprocess_function)
        self.preprocess_functions[modelName] = preprocess_function

    def Query(self, chr, start, end, modelName):
        pred_vals = []

        for index ,row in self.expectedMeasurement[0].get_data(chr, start, end)[0].iterrows():
            print(index)
            pred_val = self.models[modelName].predict(row.chr, row.start, row.end, {})
            pred_vals.append(np.ravel(pred_val)[0])
        # query measurement for expected
        result = {
                    
                    'model': modelName,
                    'value': np.array(pred_vals),
                    'expected_value': np.array(self.expectedMeasurement[0].get_data(chr, start, end)[0][self.expectedMeasurement[0].mid]),
                    'type': 'ml_model'
        }
        return result
    
    def compute(self, chr, start, end, params=None):
        # run predict on all models
        modelNames = self.models.keys()
        results = []
        for name in modelNames:
            result = self.Query(chr, start, end, name)
            corr, pvalue = pearsonr(result["expected_value"], result["value"])
            
            if pvalue <= self.pval_threshold:
                results.append(
                    {
                        'measurements': (self.expectedMeasurement[0], self.expectedMeasurement[0]),
                        'test': 'ml models',
                        'value': corr,
                        'pvalue': pvalue,
                        'type': 'ml_model'
                    }
                )
            # run pearson on expected value
        return self.toDataFrame(results)
