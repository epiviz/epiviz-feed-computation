import pandas as pd
from .stats import Correlation, CorrelationGeneSignal, OverlapBlock, TtestBlock, TtestExp, LondFDR
from .ml_models import ModelManager
from .ml_models import helper_functions

def computational_request(chr, start, end, gene_name, measurements, computations=None, model_params=None, pval_threshold=None):
    if computations is None:
        computations = ["Correlation", "CorrelationGeneSignal", "OverlapBlock", "TtestBlock", "TtestExp", "ModelManager"]
    
    for comp in computations:
        if comp != "ModelManager":
            compObj = eval(comp)(measurements, pval_threshold)
            yield compObj.compute(chr, start, end, params={"datatype": "expression", "annotation":None})

    modelObj = eval("ModelManager")(measurements, pval_threshold)
    for model_param in model_params:
        getDataFunc = getattr(helper_functions, model_param[2])
        modelObj.Add(model_param[0], model_param[1], getDataFunc, model_param[3])

    yield modelObj.compute(chr, start, end)