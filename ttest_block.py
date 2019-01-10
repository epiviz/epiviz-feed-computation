import json
import pandas as pd
from scipy.stats import ttest_ind
from utils import build_obj, format_expression_block_data

class Ttest_Block:

    def __init__(self, exp_data, block_data, exp_datasource, datasource_types):
        self.exp_data = exp_data
        self.block_data = block_data
        self.exp_datasource = exp_datasource
        self.datasource_types = datasource_types

        self.exp_data['index_col'] = self.exp_data.index
        self.gene_expression_block = dict()
        self.gene_expression_nonblock = dict()

    def get_expressions(self, row, exp_types, block_type):

        start = row["start"]
        end = row["end"]
        exp_srt = self.exp_data["start"]
        exp_end = self.exp_data["end"]

        exp_block = pd.DataFrame(columns=exp_types)
        in_block = ((exp_srt <= end) & (exp_end >= start)) | ((exp_end >= start) & (exp_end <= end)) | ((exp_srt >= start) & (exp_srt <= end))

        #queries the dataframe where expressions are in/overlap blocks and drops nan values
        exp_indices = list((self.exp_data.where(in_block)['index_col']).dropna().unique())

        exp_block = self.exp_data.iloc[exp_indices][exp_types]
        exp_nonblock = self.exp_data[(exp_end < start) | (exp_srt > end)][exp_types]

        self.gene_expression_block[block_type] = self.gene_expression_block.get(block_type, pd.DataFrame(columns=exp_types)).append(exp_block)
        self.gene_expression_nonblock[block_type] = self.gene_expression_nonblock.get(block_type, pd.DataFrame(columns=exp_types)).append(exp_nonblock)

    def sort_expression_blocks(self, block_type, block_dataframe):
        if not block_dataframe.empty:
            tissue_type = block_type.split("_")[1]
            exp_types = [tissue_type + "___normal", tissue_type + "___tumor"]
            block_dataframe.apply(lambda row: self.get_expressions(row, exp_types, block_type), axis=1)

    def ttest(self):
        for block_type, block_dataframe in self.block_data.items():
            self.sort_expression_blocks(block_type, block_dataframe)

        ttest_res = []
        pd_block = pd.DataFrame(self.datasource_types)
        pd_expression = pd.DataFrame(self.exp_datasource)

        for block_type, gene_per_block_exp in self.gene_expression_block.items():
            gene_per_nonblock_exp = self.gene_expression_nonblock[block_type]

            for exp_type in gene_per_block_exp:
                gene_block_exp = gene_per_block_exp[exp_type]

                if not gene_block_exp.empty:

                    gene_nonblock_exp = gene_per_nonblock_exp[exp_type]
                    t_value, p_value = ttest_ind(gene_block_exp, gene_nonblock_exp, equal_var=False)

                    print("block:" + block_type + ", gene:" + exp_type)
                    print(p_value)

                    gene_ds = json.loads(pd_expression.loc[pd_expression['id'] == exp_type].to_json(orient='records')[1: -1])
                    block_ds = json.loads(pd_block.loc[pd_block['id'] == block_type].to_json(orient='records')[1: -1])

                    data = format_expression_block_data(gene_block_exp, gene_nonblock_exp)
                    ttest_obj = build_obj('t-test', 'expression', 'block', False, gene_ds, block_ds, t_value, p_value, data)

                    ttest_res.append(ttest_obj)

        ttest_res = sorted(ttest_res, key=lambda x: x['value'], reverse=True)
        return ttest_res
