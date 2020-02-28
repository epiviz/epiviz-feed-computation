# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound
from .stats import Correlation, CorrelationGeneSignal, OverlapBlock, TtestBlock, TtestExp, LondFDR
from .interface import computational_request
try:
    # Change here if project is renamed and does not equal the package name
    dist_name = 'epivizFeedCompute'
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound
