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
        self.expectedMeasurements = {}
        self.pval_threshold = pval_threshold
        self.measurements = measurements

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

    def Add(self, modelName, model, preprocess_function, expected_type):
        self.models[modelName] = BaseModel(self.measurements, model, preprocess_function)
        self.preprocess_functions[modelName] = preprocess_function
        self.expectedMeasurements[modelName] = self.query_measurements(self.measurements, expected_type)

    def Query(self, chr, start, end, modelName):
        # it returns an 1x1 2D array
        pred_val = self.models[modelName].predict(chr, start, end, {})[0][0]
        exp_meas = self.expectedMeasurements[modelName][0]
        exp_val = np.mean(np.array(exp_meas.get_data(chr, start, end)[0][exp_meas.mid]))
        result = {
                    
                    'model': modelName,
                    'value': pred_val,
                    'expected_value': exp_val,
                    'type': 'ml_model'
        }
        return result
    
    def compute(self, chr, start, end, params=None):
        # run predict on all models
        modelNames = self.models.keys()
        results = []
        output = pd.DataFrame()
        for name in modelNames:
            output = output.append(self.Query(chr, start, end, name),ignore_index=True)
        corr, pvalue = pearsonr(output["expected_value"], output["value"])
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
