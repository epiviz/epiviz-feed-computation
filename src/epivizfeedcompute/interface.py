import pandas as pd
from .stats import Correlation, CorrelationGeneSignal, OverlapBlock, TtestBlock, TtestExp, LondFDR
from .ml_models import ModelManager

def computational_request(chr, start, end, gene_name, measurements, computations=None, pval_threshold=None):
    if computations is None:
        computations = ["Correlation", "CorrelationGeneSignal", "OverlapBlock", "TtestBlock", "TtestExp", "ModelManager"]
                
    for comp in computations:
        compObj = eval(comp)(measurements, pval_threshold)
        yield compObj.compute(chr, start, end, params={"datatype": "expression", "annotation":None})
