import json
import math
import pandas as pd
from scipy.stats import ttest_ind, norm
from requests import get_sample_counts
from utils import build_obj, format_expression_block_data, build_exp_singlegene_obj


class Ttest_Gene:

    def __init__(self, gene_types, exp_data, chromosome, start_seq, end_seq):
        self.gene_types = gene_types
        self.exp_data = exp_data
        self.chromosome = chromosome
        self.start_seq = start_seq
        self.end_seq = end_seq

    def ttest_calculations(self, row, gene1, gene2, sample_counts):
        gene_norm = gene1["id"]
        gene_cancer = gene2["id"]

        one = row[gene_norm]
        two = row[gene_cancer]

        variance_threshold = 0.05 * 0.95

        var_one = variance_threshold if (one * (1 - one)) < variance_threshold else (one * (1 - one))

        var_two = variance_threshold if (two * (1 - two)) < variance_threshold else (two * (1 - two))

        denominator = math.sqrt(var_one / sample_counts[gene_norm][0] + var_two / sample_counts[gene_cancer][0])

        ttest_value = (one - two) / denominator

        p_value = 1 - norm.cdf(ttest_value)

        data = [{
            "type": gene1["name"],
            "value": one
        }, {
            "type": gene2["name"],
            "value": two
        }]

        corr_obj = [build_exp_singlegene_obj('Binomial test difference in proportions', 'expression', 'expression',
                    True, gene1, gene2, ttest_value, pvalue=p_value, gene=row['gene'], data=data)]
        return corr_obj

    def ttest_exp(self):
        print ("ttest per single gene!")
        sample_counts = get_sample_counts(self.gene_types, self.start_seq, self.end_seq, self.chromosome)

        ttest_results = []
        if self.exp_data.empty or not sample_counts:
            return ttest_results
        # sorts gene_types by tissue type
        sorted_gene_types = sorted(self.gene_types, key=lambda x: x["id"])
        # gets list of dicts containing the same info as gene type
        tissue_types = [item for item in sorted_gene_types]
        # breaks list down into sublist of size 2
        gene_pairs = [tissue_types[x:x+2] for x in range(0, len(tissue_types), 2)]

        for gene_pair in gene_pairs:
            # preforms ttest calc on each row of the exp_data dataframe
            results = self.exp_data.apply(lambda row: self.ttest_calculations(row, gene_pair[0], gene_pair[1], sample_counts), axis=1)
            # concat series of list
            concat_result = [element for sublist in results for element in sublist]
            ttest_results.extend(concat_result)

        ttest_results = sorted(ttest_results, key=lambda x: x['value'], reverse=True)
        return ttest_results
