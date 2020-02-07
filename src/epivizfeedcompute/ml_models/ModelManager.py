import tensorflow as tf
import keras
import helper_functions
import pandas as pd
import numpy as np
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Activation, Flatten, Conv1D, MaxPooling1D
from keras.layers.core import Reshape
from .BaseModel import BaseModel

class ModelManager(BaseModel):

    def __init__ (self, measurements):
        self.models = {}
        self.measurements = measurements

    def Add(self, modelName, model):
        self.models[modelName] = model

    def Query(chr, start, end, modelName):
        return self.models[modelName].predict(chr, start, end)