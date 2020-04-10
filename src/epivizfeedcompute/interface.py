import pandas as pd
from .stats import Correlation, CorrelationGeneSignal, OverlapBlock, TtestBlock, TtestExp, LondFDR
from .ml_models import ModelManager
from .ml_models import helper_functions

def computational_request(chr, start, end, gene_name, measurements, computations=None, models=None, pval_threshold=None):
    if computations is None:
        computations = ["Correlation", "CorrelationGeneSignal", "OverlapBlock", "TtestBlock", "TtestExp", "ModelManager"]
    
    for comp in computations:
        compObj = eval(comp)(measurements, pval_threshold)
        yield compObj.compute(chr, start, end, params={"datatype": "expression", "annotation":None})

    model_params = ["model1", "epiviz_model_1.model", helper_functions.get_data]
    modelObj = eval("ModelManager")(measurements, pval_threshold)
    for model in models:
        modelObj.Add(model_params[0], model_params[1], model_params[2])

    yield modelObj.compute(chr, start, end)