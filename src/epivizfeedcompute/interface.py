import pandas as pd
from .stats import Correlation, CorrelationGeneSignal, OverlapBlock, TtestBlock, TtestExp, LondFDR
from .ml_models import ModelManager
from .ml_models import helper_functions
def computational_request(chr, start, end, gene_name, measurements, computations=None, model_params = None,pval_threshold=None):
    if computations is None:
        computations = ["Correlation", "CorrelationGeneSignal", "OverlapBlock", "TtestBlock", "TtestExp", "ModelManager"]
    if model_params is None:
        model_params = ["model1","..\ml_models\epiviz_model_1.model",helper_functions.get_data]
    for comp in computations:
        compObj = eval(comp)(measurements, pval_threshold)
        if comp == "ModelManager":
            compObj.Add(model_params[0], model_params[1], model_params[2])
        yield compObj.compute(chr, start, end, params={"datatype": "expression", "annotation":None})
