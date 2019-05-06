import numpy as np
import pandas as pd
import logging
import math
import json
import itertools
from old_feed.utils import build_obj
from scipy.stats import ttest_ind, fisher_exact, norm
from epivizFeed.StatMethod import StatMethod
from old_feed.data_functions import Block_data


class OverlapBlock(StatMethod):
    def __init__(self, measurements):
        super(OverlapBlock, self).__init__(measurements)
        self.data_sources = super(OverlapBlock, self).get_measurements_self("block")

    def calc_region(self, block, start_seq, end_seq):
        for start, end in zip(block["block_tissue"]['start'], block["block_tissue"]['end']):
            min_end = min(end, float(end_seq))
            max_srt = max(start, float(start_seq))
            if min_end > max_srt:
                block['region'].append(min_end - max_srt)

    def get_overlap(self, overlap_region, block_one, block_two, start_seq, end_seq):
        b_one_ind = block_one["index"]
        b_two_ind = block_two["index"]

        while b_one_ind < block_one["length"] and b_two_ind < block_two["length"]:

            tissue_one_start = max(float(start_seq), block_one["block_tissue"]['start'][b_one_ind])
            tissue_two_start = max(float(start_seq), block_two["block_tissue"]['start'][b_two_ind])
            tissue_one_end = min(float(end_seq), block_one["block_tissue"]['end'][b_one_ind])
            tissue_two_end = min(float(end_seq), block_two["block_tissue"]['end'][b_two_ind])
            # there is an overlap
            if tissue_one_start <= tissue_two_start < tissue_one_end or \
               tissue_one_start < tissue_two_end <= tissue_one_end or \
               (tissue_one_start == tissue_two_start and tissue_one_end == tissue_two_end):

                common_end = min(tissue_two_end, tissue_one_end)
                common_start = max(tissue_one_start, tissue_two_start)
                if common_end > common_start:
                    overlap_region.append(common_end - common_start)
            # block tissue two is larger
            if tissue_two_start >= tissue_one_end or tissue_two_end > tissue_one_end:
                b_one_ind += 1
            # block tissue one is larger
            if tissue_one_start >= tissue_two_end or tissue_two_end <= tissue_one_end:
                b_two_ind += 1

    def calc_overlap_percentage(self, overlap_region, block_one, block_two, data_source_one, data_source_two, start_seq, end_seq):
        self.calc_region(block_one, start_seq, end_seq)
        self.calc_region(block_two, start_seq, end_seq)
        self.get_overlap(overlap_region, block_one, block_two, start_seq, end_seq)

        sum_block_one_region = sum(block_one["region"])
        sum_block_two_region = sum(block_two["region"])
        overlap = sum(overlap_region)
        union = sum_block_one_region + sum_block_two_region - overlap
        block_one_only = max(sum_block_one_region - overlap, 0)
        block_two_only = max(sum_block_two_region - overlap, 0)
        non_block = max(int(end_seq) - int(start_seq) - union, 0)

        fisher_table = np.array([[overlap, block_one_only], [block_two_only, non_block]])
        odds_ratio, p_value = fisher_exact(fisher_table)
        if not math.isnan(odds_ratio):
            logging.info('p value is ' + str(p_value))
            logging.info('odds ratio is ' + str(odds_ratio))

            overlap_percent = 0.0 if union == 0.0 else overlap * 1.0 / union
            overlap_obj = build_obj('overlap', 'block', 'block', False, data_source_one, data_source_two, overlap_percent, p_value)

        return overlap_obj

    def create_block(self, attributes, attributes_vals):
        return {key: val for key, val in zip(attributes, attributes_vals)}

    def compute(self, chromosome, start, end, additional=None):
        block_data = Block_data(start, end, chromosome, measurements=self.data_sources)
        overlaps = []
        if block_data:
            for data_source_one, data_source_two in itertools.combinations(self.data_sources, 2):
                tissue_one = data_source_one["id"]
                tissue_two = data_source_two["id"]

                if tissue_one in block_data and tissue_two in block_data:

                    attributes = ["tissue", "block_tissue", "index", "length", "region"]
                    block_one = self.create_block(attributes, [tissue_one,
                        block_data[tissue_one], 0, len(block_data[tissue_one]['start']), []])

                    block_two = self.create_block(attributes, [tissue_two,
                        block_data[tissue_two], 0, len(block_data[tissue_two]['start']), []])

                    overlap_region = []
                    overlap_obj = self.calc_overlap_percentage(overlap_region, block_one, block_two, data_source_one, data_source_two, start, end)
                    overlaps.append(overlap_obj)
        block_overlap = pd.Series(overlaps)
        block_overlap = block_overlap.apply(pd.Series)
        block_overlap = block_overlap.to_json(orient='records')
        parse_res = json.loads(block_overlap)

        logging.info('overlap done!')
        return parse_res
