import json
import math
import pandas as pd
from scipy.stats import ttest_ind, norm
from old_feed.requests import get_sample_counts
from old_feed.utils import build_obj, format_expression_block_data, build_exp_singlegene_obj
from stat_classes.StatMethod import StatMethod
from old_feed.data_functions import Gene_data


class TtestGene(StatMethod):

    def __init__(self, measurements):
        super().__init__(measurements)
        self.measurements = measurements
        self.gene_types = super().get_measurements_self("gene")

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

        corr_obj = build_exp_singlegene_obj('Binomial test difference in proportions', 'expression', 'expression',
                    True, gene1, gene2, ttest_value, pvalue=p_value, gene=row['gene'], data=data)
        return corr_obj

    def partion(self, attr):
        #attributes other than tissues must be specifed in form [attr_name, (attr_val 1, attr_val 2)]
        # where attributes values must be distinct strings
        pair = None
        m = pd.DataFrame()
        for measurement in self.measurements:
            m = m.append(measurement, ignore_index=True)
        if attr[0] == "tissue":
            exp_measurements = m[m['name'].str.contains('Expression')]
            norm = exp_measurements[exp_measurements['name'].str.contains('normal')]
            canc = exp_measurements[exp_measurements['name'].str.contains('tumor')]
            pair = (norm, canc)
        else:
            attr_v1 = attr[1][0]
            attr_v2 = attr[1][1]

            exp_measurements = m[m['name'].str.contains('Expression')]
            attr_val_one = exp_measurements[exp_measurements['annotation'].str.contains('attr_v1')]
            attr_val_two = exp_measurements[exp_measurements['annotation'].str.contains('attr_v2')]
            pair = (attr_val_one, attr_val_two)

        return pair

    def compute(self, chromosome, start, end):
        exp_data = Gene_data(start, end, chromosome, measurements=self.gene_types)
        print("ttest per single gene!")
        sample_counts = get_sample_counts(self.gene_types, start, end, chromosome)
        cols = ['computation-type', 'data-type-one', 'data-type-two', 'show-chart', 'attribute-one', 'attribute-two', 'value', 'pvalue', 'data', 'gene']
        ttest_results = pd.DataFrame(columns=cols)

        group_one, group_two = self.partion(["tissue", (None, None)])
        # use to_dict('records') convert dataframe into a list of dicts for gene data function
        exp_group_one = Gene_data(start, end, chromosome, measurements=group_one.to_dict('records'))
        exp_group_two = Gene_data(start, end, chromosome, measurements=group_two.to_dict('records'))
        #for ele_one, ele_two in itertools.combinations()

        if exp_data.empty or not sample_counts:
            return ttest_results

        # sorts gene_types by tissue type
        sorted_gene_types = sorted(self.gene_types, key=lambda x: x["id"])
        # gets list of dicts containing the same info as gene type
        tissue_types = [item for item in sorted_gene_types]
        # breaks list down into sublist of size 2
        gene_pairs = [tissue_types[x:x+2] for x in range(0, len(tissue_types), 2)]
        print(gene_pairs)
        for gene_pair in gene_pairs:
            # preforms ttest calc on each row of the exp_data dataframe
            results = exp_data.apply(lambda row: self.ttest_calculations(row, gene_pair[0], gene_pair[1], sample_counts), axis=1)
            # turns series of dicts to a dataframe
            results = results.apply(pd.Series)

            ttest_results = pd.concat([ttest_results, results])
        ttest_results = ttest_results.sort_values(by=['value'])
        # print(ttest_results["value"])
        print("ttest_gene_result")
        return ttest_results
