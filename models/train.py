import tensorflow as tf
import pickle
import keras
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from sklearn.linear_model import LinearRegression
from epivizfileserver.measurements import WebServerMeasurement
from epivizfeedcompute.stats import BaseStats, TtestBlock, TtestExp, Correlation, CorrelationGeneSignal, OverlapBlock
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv1D, MaxPooling1D
from keras.layers.core import Reshape
from sklearn.model_selection import train_test_split


dataF = pickle.load( open( "expression_lung_data.p", "rb" ) )
in_gene = pickle.load( open( "in_gene_lung_tumor_buckets_v2.p", "rb" ) )
in_gene_islands = pickle.load( open( "in_gene_lung_cpg_islands.p", "rb" ) )
upstream_islands = pickle.load( open( "methy_upstream_lung_cpg_islands.p", "rb" ) )
downstream_islands = pickle.load( open( "methy_downstream_lung_cpg_islands.p", "rb" ) )
methy_upstream = pickle.load( open( "methy_upstream_lung_tumor_buckets_v2.p", "rb" ) )
methy_downstream = pickle.load( open( "methy_downstream_lung_tumor_buckets_v2.p", "rb" ) )

col_names = []
for i in range(20):
    col_names.append("in_gene {}".format(i))
in_gene_df = pd.DataFrame(in_gene, columns = col_names)

col_names = []
for i in range(10):
    col_names.append("upstream {}".format(i))
upstream_df = pd.DataFrame(methy_upstream, columns = col_names)

col_names = []
for i in range(10):
    col_names.append("downstream {}".format(i))
downstream_df = pd.DataFrame(methy_downstream, columns = col_names)

cpg_islands_in_gene = pd.DataFrame(in_gene_islands,columns = ["cpg_count", "contains_cpg_islands"])
cpg_islands_upstream = pd.DataFrame(upstream_islands,columns = ["upstream_cpg_count"])
cpg_islands_downstream = pd.DataFrame(downstream_islands,columns = ["downstream_cpg_count"])

in_gene_df = in_gene_df.fillna(0)
upstream_df = upstream_df.fillna(0)
downstream_df = downstream_df.fillna(0)
X = pd.concat([in_gene_df,upstream_df,downstream_df,cpg_islands_in_gene,cpg_islands_upstream,cpg_islands_downstream], axis=1)
y = pd.Series(dataF["lung___tumor"])

X_train, X_test, y_train, y_test = train_test_split( X, y, test_size=0.2)
n_timesteps, n_features, n_outputs = X_train.shape[0], X_train.shape[1], y_train.shape[0]
verbose, epochs, batch_size = 0, 100, 85

verbose, epochs, batch_size = 0, 100, 85

model = Sequential()
model.add(Conv1D(filters=64, kernel_size=10, activation='relu', input_shape=(n_features, 1)))
model.add(MaxPooling1D(pool_size=3))
model.add(Flatten())
model.add(Dense(600, activation='relu'))
model.add(Dropout(0.1))
model.add(Dense(250, activation='relu'))
# model.add(Dense(125, activation='relu'))
model.add(Dense(1, activation='relu'))
model.compile(loss='mean_squared_error', optimizer='sgd', metrics=['mean_squared_error'])
print(model.summary())
X_train = np.expand_dims(X_train, axis=2)
X_test = np.expand_dims(X_test, axis=2)
y_train = np.array(y_train)
# X_train = np.reshape(X_train,(1, X_train.shape[0], X_train.shape[1]))
print(X_train.shape)
model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, verbose=verbose)
_, accuracy = model.evaluate(X_test, y_test, batch_size=batch_size, verbose=0)

model.save('epiviz_model_lung.model')