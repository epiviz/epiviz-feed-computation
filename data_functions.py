import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind, fisher_exact, norm
from scipy.io import savemat
from utils import build_obj, build_exp_methy_obj, add_to_block, format_expression_block_data, build_exp_singlegene_obj
from requests import get_methy_data, get_block_data, get_gene_data, get_sample_counts
from urllib.request import urlopen
import json
import itertools
import math

class Data_Functions:

    def __init__(self, start_seq, end_seq, chromosome, gene_name, measurements = None):
        self.start_seq = start_seq
        self.end_seq = end_seq
        self.chromosome = chromosome
        self.gene_name = gene_name
        self.measurements = measurements

    def gene_data(self):
        return get_gene_data(self.start_seq, self.end_seq, self.chromosome, self.measurements)

    def block_data(self):

        return get_block_data(self.start_seq, self.end_seq, self.chromosome, self.measurements)

    def methylation_diff(self):
        return get_methy_data(self.start_seq, self.end_seq, self.chromosome, self.measurements)

    def methylation(self):
        return get_methy_data(self.start_seq, self.end_seq, self.chromosome, self.measurements)
