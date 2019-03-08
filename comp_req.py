import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind, fisher_exact, norm
from scipy.io import savemat
from old_feed.utils import build_obj, build_exp_methy_obj, add_to_block, format_expression_block_data, build_exp_singlegene_obj
from old_feed.requests import get_methy_data, get_block_data, get_gene_data, get_sample_counts
from old_feed.data_functions import Gene_data, Block_data, Methylation_diff, Methylation
from urllib.request import urlopen
from stat_classes.ttest_block import TtestBlock
from stat_classes.ttest_gene import TtestGene
from stat_classes.overlap_block_percent import OverlapBlock
from stat_classes.correlation_exp_methy import CorrelationExpMethy
from stat_classes.correlation_exp import CorrelationExp
from stat_classes.correlation_methy import CorrelationMethy
import json
import itertools
import math


def comp_req(start_seq, end_seq, chromosome, gene_name, measurements=None):

    ttest_gene = TtestGene(measurements)
    yield ttest_gene.compute(chromosome, start_seq, end_seq)

    block_ol = OverlapBlock(measurements)
    yield block_ol.compute(chromosome, start_seq, end_seq)

    methy_diff = CorrelationMethy(measurements, "methy_diff")
    yield methy_diff.compute(chromosome, start_seq, end_seq)

    methy_corr = CorrelationMethy(measurements, "methy")
    yield methy_corr.compute(chromosome, start_seq, end_seq)

    corr_exp = CorrelationExp(measurements)
    yield corr_exp.compute(chromosome, start_seq, end_seq)

    ttest = TtestBlock(measurements)
    yield ttest.compute(chromosome, start_seq, end_seq)

    test_method = CorrelationExpMethy(measurements, "methy")
    yield test_method.compute(chromosome, start_seq, end_seq)

    test_method = CorrelationExpMethy(measurements, "methy_diff")
    yield test_method.compute(chromosome, start_seq, end_seq)
