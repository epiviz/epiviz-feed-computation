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
