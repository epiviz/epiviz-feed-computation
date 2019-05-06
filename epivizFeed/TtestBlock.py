import json
import pandas as pd
import logging
from scipy.stats import ttest_ind
from old_feed.utils import build_obj, format_expression_block_data
from epivizFeed.StatMethod import StatMethod
from old_feed.data_functions import Gene_data, Block_data


class TtestBlock(StatMethod):

    def __init__(self, measurements):
        super(TtestBlock, self).__init__(measurements)
        self.exp_datasource = super(TtestBlock, self).get_measurements_self("gene")
        self.datasource_types = super(TtestBlock, self).get_measurements_self("block")

    def get_expressions(self, row, exp_types, block_type, exp_data, exp_block, exp_nonblock):
        start = row["start"]
        end = row["end"]
        exp_srt = exp_data["start"]
        exp_end = exp_data["end"]

        blocks = pd.DataFrame(columns=exp_types)
        # boolean formula for finding expression block overlap
        in_block = ((exp_srt <= end) & (exp_end >= start)) | ((exp_end >= start) & (exp_end <= end)) | ((exp_srt >= start) & (exp_srt <= end))
        # queries the dataframe where expressions are in/overlap blocks and drops nan values
        exp_indices = list((exp_data.where(in_block)['index_col']).dropna().unique())
        # gets rows at exp_indices keeping only the exp types cols
        blocks = exp_data.iloc[exp_indices][exp_types]

        nonblocks = exp_data[(exp_end < start) | (exp_srt > end)][exp_types]

        exp_block[block_type] = exp_block.get(block_type, pd.DataFrame(columns=exp_types)).append(blocks)
        exp_nonblock[block_type] = exp_nonblock.get(block_type, pd.DataFrame(columns=exp_types)).append(nonblocks)

    def ttest_calculation(self, gene_block_exp, gene_per_nonblock_exp, exp_type, block_type, pd_block, pd_expression):

        gene_nonblock_exp = gene_per_nonblock_exp[exp_type]
        t_value, p_value = ttest_ind(gene_block_exp, gene_nonblock_exp, equal_var=False)

        print("block:" + block_type + ", gene:" + exp_type)
        print(p_value)

        gene_ds = json.loads(pd_expression.loc[pd_expression['id'] == exp_type].to_json(orient='records')[1: -1])
        block_ds = json.loads(pd_block.loc[pd_block['id'] == block_type].to_json(orient='records')[1: -1])

        data = format_expression_block_data(gene_block_exp, gene_nonblock_exp)
        ttest_obj = build_obj('t-test', 'expression', 'block', False, gene_ds, block_ds, t_value, p_value, data)

        return ttest_obj

    def partition_data(self, block_data, exp_data, add_attr=None):
        exp_block = dict()
        exp_nonblock = dict()

        for block_type, block_dataframe in block_data.items():
            if not block_dataframe.empty:
                tissue_type = block_type.split("_")[1]
                exp_types = [tissue_type + "___normal", tissue_type + "___tumor"]
                # gets block and non-block expressions
                block_dataframe.apply(lambda row: self.get_expressions(row, exp_types, block_type, exp_data, exp_block, exp_nonblock), axis=1)

        return 1

    def compute(self, chromosome, start, end, additional=None):
        exp_data = Gene_data(start, end, chromosome, measurements=self.exp_datasource)
        block_data = Block_data(start, end, chromosome, measurements=self.datasource_types)
        exp_data['index_col'] = exp_data.index
        gene_expression_block = dict()
        gene_expression_nonblock = dict()
        # loop through block of different tissue types
        for block_type, block_dataframe in block_data.items():
            if not block_dataframe.empty:
                tissue_type = block_type.split("_")[1]
                exp_types = [tissue_type + "___normal", tissue_type + "___tumor"]
                # gets block and non-block expressions
                block_dataframe.apply(lambda row: self.get_expressions(row, exp_types, block_type, exp_data, gene_expression_block, gene_expression_nonblock), axis=1)

        ttest_res = pd.DataFrame()
        results = []
        pd_block = pd.DataFrame(self.datasource_types)
        pd_expression = pd.DataFrame(self.exp_datasource)

        for block_type, gene_per_block_exp in gene_expression_block.items():
            exp_types = list(gene_per_block_exp.columns)
            gene_per_nonblock_exp = gene_expression_nonblock[block_type]

            for exp_type in exp_types:
                gene_block_exp = gene_per_block_exp[exp_type]

                if not gene_block_exp.empty:
                    gene_nonblock_exp = gene_per_nonblock_exp[exp_type]
                    ttest_obj = self.ttest_calculation(gene_block_exp, gene_per_nonblock_exp, exp_type, block_type, pd_block, pd_expression)
                    results.append(ttest_obj)

        results = sorted(results, key=lambda x: x['value'], reverse=True)
        ttest_res = pd.Series(results)
        ttest_res = ttest_res.apply(pd.Series)
        ttest_res = ttest_res.to_json(orient='records')
        parse_res = json.loads(ttest_res)

        logging.info("ttest_block_res")
        return parse_res
