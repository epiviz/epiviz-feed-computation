import tensorflow as tf
import keras
import pandas as pd
import numpy as np
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Activation, Flatten, Conv1D, MaxPooling1D
from keras.layers.core import Reshape
from .BaseModel import BaseModel

class ModelManager():

    def __init__ (self, measurements, pval_threshold):
        self.models = {}
        self.preprocess_functions = {}
        self.measurements = measurements

    def Add(self, modelName, model, preprocess_function):
        self.models[modelName] = BaseModel(self.measurements, model, preprocess_function)
        self.preprocess_functions[modelName] = preprocess_function

    def Query(self, chr, start, end, modelName):
        pred_vals = self.models[modelName].predict(chr, start, end, {})
        print("dsfa")
        print(pred_vals)
        result = {
                            'model': modelName,
                            'value': pred_vals,
                            'type': 'ml_model'
        }
        return result
    
    def compute(self, chr, start, end, params=None):
        # run predict on all models
        modelNames = self.models.keys()
        results = []
        for name in modelNames:
            results.append(self.Query(chr, start, end, name))
        return pd.DataFrame(results)

        #     for (m1, m2) in msets:
        #     data1, data2 = self.get_transform_data([m1, m2], chr, start, end, params)
        #     corr, pvalue = self.compute_stat(data1, data2)
            
        #     if pvalue <= self.pval_threshold:
        #         results.append(
        #             {
        #                 'measurements': (m1, m2),
        #                 'test': 'correlation',
        #                 'value': corr,
        #                 'pvalue': pvalue,
        #                 'type': 'computation'
        #             }
        #         )

        # sorted_results = sorted(results, key=lambda x: x['value'], reverse=True)
        # return pd.DataFrame(
        #     data = {
        #         "measurement1" : measurement1,
        #         "measurement2" : measurement2,
        #         "value" : value,
        #         "pvalue" : pvalue
        #     }
        # )

