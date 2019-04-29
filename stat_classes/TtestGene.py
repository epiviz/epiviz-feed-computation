import json
import math
import itertools
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

    def partion(self, type, group_one, group_two=None):
        #attributes other than tissues must be specifed in form [attr_name, (attr_val 1, attr_val 2)]
        # where attributes values must be distinct strings
        pair = None
        m = pd.DataFrame()
        for measurement in self.measurements:
            m = m.append(measurement, ignore_index=True)
        exp_measurements = m[m['name'].str.contains('Expression')]
        print(m)

        if group_two is None:
            g_one = exp_measurements[exp_measurements['name'].str.contains(group_one)]
            g_two = exp_measurements[exp_measurements['name'].str.contains(group_one) == False]
            pair = (g_one, g_two)
        else:
            g_one = exp_measurements[exp_measurements['annotation'].str.contains(group_one)]
            g_two = exp_measurements[exp_measurements['annotation'].str.contains(group_two)]
            pair = (g_one, g_two)
        return pair

    def to_list_of_dict(self, ele):
        ret_val = None
        for i in self.gene_types:
            if i["id"] == ele:
                ret_val = i

        return ret_val

    def grouping(self, group_one, group_two, all_pairs=True):
        group_pairs = []

        if all_pairs:

            g_one = [c for c in group_one.columns if "_" in c]
            g_two = [c for c in group_two.columns if "_" in c]
            group_pairs = [(x, y) for x in g_one for y in g_two]
            group_pairs = [(self.to_list_of_dict(x), self.to_list_of_dict(y)) for x, y in group_pairs]
        else:
            a = group_one.columns
            b = group_two.columns

            group_pairs = [(x, y) for x, y in zip(a, b) if x not in b and y not in a]
            group_pairs = [(self.to_list_of_dict(x), self.to_list_of_dict(y)) for x, y in group_pairs]

        return group_pairs

    def unpack_params(self, additional):
        if additional is not None:
            partition_type = additional["partition_type"]
            group_one = additional["group_one"]
            group_two = additional["group_two"]
            grouping = additional["grouping"]
        else:
            partition_type = "condition"
            group_one = "normal"
            group_two = "tumor"
            grouping = "one-to-one"
        return partition_type, group_one, group_two, grouping

    def compute(self, chromosome, start, end, additional=None):
        part_type, g_one, g_two, grouping = self.unpack_params(additional)
        if grouping == "all_pairs":
            grouping = True
        else:
            grouping = False

        exp_data = Gene_data(start, end, chromosome, measurements=self.gene_types)
        print("ttest per single gene!")
        group_one, group_two = self.partion(part_type, g_one)
        exp_group_one = Gene_data(start, end, chromosome, measurements=group_one.to_dict('records'))
        exp_group_two = Gene_data(start, end, chromosome, measurements=group_two.to_dict('records'))

        group_pairs = self.grouping(exp_group_one, exp_group_two, grouping)
        print(group_pairs)
        sample_counts = get_sample_counts(self.gene_types, start, end, chromosome)
        cols = ['computation-type', 'data-type-one', 'data-type-two', 'show-chart', 'attribute-one', 'attribute-two', 'value', 'pvalue', 'data', 'gene']
        ttest_results = pd.DataFrame(columns=cols)
        if exp_data.empty or not sample_counts:
            return ttest_results

        for g_pair in group_pairs:

            # preforms ttest calc on each row of the exp_data dataframe
            results = exp_data.apply(lambda row: self.ttest_calculations(row, g_pair[0], g_pair[1], sample_counts), axis=1)
            # turns series of dicts to a dataframe
            results = results.apply(pd.Series)

            ttest_results = pd.concat([ttest_results, results])

        ttest_results = ttest_results.sort_values(by=['value'], ascending=False)
        ttest_results = ttest_results[cols]
        print(ttest_results)
        # front end assumes that json parse returns string but when takes in a json
        # string ret a str so we turn ttest_result to a json string then to a list for proper behavior
        # turns data into a json string
        ttest_results = ttest_results.to_json(orient='records')
        parse_res = json.loads(ttest_results)

        print("ttest_gene_result")

        return parse_res
