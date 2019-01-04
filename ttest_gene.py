import json
import pandas as pd
from scipy.stats import ttest_ind
from utils import build_obj, format_expression_block_data


class ttest_block:

    def __init__(self, gene_types, exp_data, chromosome, start_seq, end_seq):
        self.gene_types = gene_types
        self.exp_data = exp_data
        self.chromosome = chromosome
        self.start_seq = start_seq
        self.end_seq = end_seq

    def ttest_expression_per_gene(gene_types, exp_data, chromosome, start_seq, end_seq):
        print ("ttest per single gene!")

        sample_counts = get_sample_counts(
            gene_types, start_seq, end_seq, chromosome)

        ttest_results = []

        if exp_data.empty or not sample_counts:
            return ttest_results

        gene_pairs = [["breast___normal", "breast___tumor"],
                      ['colon___normal', 'colon___tumor'],
                      ['lung___normal', 'lung___tumor'],
                      ['thyroid___normal', 'thyroid___tumor']]
        for gene_pair in gene_pairs:
            exp1 = gene_pair[0]
            exp2 = gene_pair[1]
            data_source_one = [
                element for element in gene_types if element["id"] == exp1][0]
            data_source_two = [
                element for element in gene_types if element["id"] == exp2][0]

            for index, row in exp_data.iterrows():

                one = row[exp1]
                two = row[exp2]

                variance_threshold = 0.05 * 0.95

                var_one = variance_threshold if (
                    one * (1 - one)) < variance_threshold else (one * (1 - one))

                var_two = variance_threshold if (
                    two * (1 - two)) < variance_threshold else (two * (1 - two))

                denominator = math.sqrt(var_one / sample_counts[exp1][0] +
                                        var_two / sample_counts[exp2][0])

                ttest_value = (one - two) / denominator

                p_value = 1 - norm.cdf(ttest_value)

                data = [{
                    "type": data_source_one["name"],
                    "value": one
                }, {
                    "type": data_source_two["name"],
                    "value": two
                }]

                corr_obj = build_exp_singlegene_obj('Binomial test difference in proportions', 'expression',                                        'expression', True, data_source_one,
                                                    data_source_two, ttest_value, pvalue=p_value, gene=row['gene'], data=data)
                ttest_results.append(corr_obj)

        ttest_results = sorted(ttest_results, key=lambda x: x['value'],
                               reverse=True)

        return ttest_results
