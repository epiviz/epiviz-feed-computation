import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind, fisher_exact, norm
from scipy.io import savemat
from old_feed.utils import build_obj, build_exp_methy_obj, add_to_block, format_expression_block_data, build_exp_singlegene_obj
from old_feed.requests import get_methy_data, get_block_data, get_gene_data, get_sample_counts

from urllib import urlopen
import json
import itertools
import math


def format_exp_methy_output(attr1, attr2, type1, type2):
    data = []
    for exp, methy in zip(attr1, attr2):
        point = dict()
        point[type1] = exp
        point[type2] = methy
        data.append(point)

    return data
