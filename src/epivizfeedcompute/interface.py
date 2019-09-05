import pandas as pd
from epivizfeedcompute.stats import CorrelationSignal


def computational_request(chr, start, end, gene_name, measurements, computations=None, pval_threshold=None):
    if computations is None:
        computations = ["CorrelationSignal"]
                
    for comp in computations:
        compObj = eval(comp)(measurements, pval_threshold)
        yield compObj.compute(chr, start, end)
