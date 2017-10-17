# This file computes correlations between different gene expression tracks.
# It was originally created to read from csv file, not used any more.
# Created in Sep. 2017 by Zhe Cui
#

import csv
import numpy as np
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind

# skip pancreas for now, since we do not have that data
tissue_types = ['breast', 'colon', 'lung', 'thyroid']

# get block data from different files

block_data = dict()
for tissue_type in tissue_types:
    block_data[tissue_type] = np.genfromtxt(
        '../data/' + tissue_type + '_blocks.csv',
        names=True,  # If `names` is True,
        # the field names are read from the first valid line
        delimiter=',',  # tab separated values
        usecols=(2, 3),
        unpack=True,
        dtype=None)  # guess the data type of each column

methy_file = '../data/methy.csv'

methy_raw = np.genfromtxt(
    methy_file,
    names=True,  # If `names` is True,
    # the field names are read from the first valid line
    delimiter=',',  # tab separated values
    usecols=(2, 3, 5, 6, 7, 8), # read in start/end ind and all methy data
    unpack=True,
    dtype=None)  # guess the data type of each column


def gene_in_block(block_tissue, gene_start, gene_end):
    for block_start, block_end in block_data[block_tissue]:
        if block_start <= gene_end <= block_end:
            return True
        if block_start <= gene_start <= block_end:
            return True
        if gene_start <= block_start <= gene_end:
            return True
        if gene_start <= block_end <= gene_end:
            return True

    return False


def gene_in_methy(gene_start, gene_end):
    for methy_start, methy_end in zip(methy_raw['start'], methy_raw['end']):
        if methy_start <= gene_end <= methy_end:
            return True
        if methy_start <= gene_start <= methy_end:
            return True
        if gene_start <= methy_start <= gene_end:
            return True
        if gene_start <= methy_end <= gene_end:
            return True

    return False


def add_to_list(block, expression, tissue_key, block_key, ind):
    key = tissue_key + '|' + block_key
    if key not in block:
        block[key] = []
    block[key].append(expression[tissue_key][ind])


def add_to_block(gene_block, gene_expression, block_type, ind):
    for tissue in tissue_types:
        add_to_list(gene_block, gene_expression,
                    tissue + '___tumor', block_type, ind)
        if tissue != 'pancreas':
            add_to_list(gene_block, gene_expression, tissue +
                        '___normal', block_type, ind)


expression_file = '../data/gene_expression_chr11.csv'
correlation_file = '../data/gene_expression_chr11_correlation.csv'

# csvfile = open(expression_file, 'r')
# csvreader = csv.reader(csvfile, delimiter=',')
write_file = open(correlation_file, 'wb')
csv_writer = csv.writer(write_file, delimiter=',',
                        quoting=csv.QUOTE_NONE)
expression_data = np.genfromtxt(
    expression_file,
    names=True,  # If `names` is True,
    # the field names are read from the first valid line
    delimiter=',',  # tab separated values
    usecols=(4, 5, 6, 7, 8, 9, 10, 11),
    unpack=True,
    dtype=None)  # guess the data type of each column

# print(type(expression_data.dtype.names))
csv_writer.writerow(['tissue type 1', 'tissue type 2', 'correlation'])
# iterate through the names of each columns
# for col1 in expression_data.dtype.names:
for ind_one in range(0, len(expression_data.dtype.names)):
    col_one = expression_data[expression_data.dtype.names[ind_one]]
    for ind_two in range(ind_one + 1, len(expression_data.dtype.names)):
        col_two = expression_data[expression_data.dtype.names[ind_two]]
        correlation_coefficient = pearsonr(col_one, col_two)
        # print correlation_coefficient
        csv_writer.writerow([expression_data.dtype.names[ind_one],
                             expression_data.dtype.names[ind_two],
                             correlation_coefficient[0]])

write_file.close()

# gene expression and block independency test
block_indi_file = open('../data/gene_block_indicator.csv', 'wb')
block_writer = csv.writer(block_indi_file, delimiter=',',
                          quoting=csv.QUOTE_NONE)
# read all the data
expression_data = np.genfromtxt(
    expression_file,
    names=True,  # If `names` is True,
    # the field names are read from the first valid line
    delimiter=',',  # tab separated values
    unpack=True,
    dtype=None)  # guess the data type of each column

# construct block indicator
block_indicator = dict()
gene_expression_block = dict()
gene_expression_nonblock = dict()
for tissue_type in tissue_types:
    block_indicator[tissue_type] = []

    gene_start_seq = expression_data['start']
    gene_end_seq = expression_data['end']
    # loop through gene expression
    for gene_ind in range(0, gene_start_seq.shape[0]):
        if gene_in_block(tissue_type, gene_start_seq[gene_ind], gene_end_seq[
           gene_ind]):
            block_indicator[tissue_type].append(1)
            add_to_block(gene_expression_block, expression_data, tissue_type,
                         gene_ind)
        else:
            block_indicator[tissue_type].append(0)
            add_to_block(gene_expression_nonblock, expression_data, tissue_type,
                         gene_ind)


# write block indicator into files, 0: false, 1: true
for index in range(0, len(block_indicator['breast'])):
    row = [block_indicator[tissue_type][index] for tissue_type in tissue_types]
    block_writer.writerow(row)

# calculate t test between block and non-block gene expression of the same
# tissue type
ttest_file = open('../data/gene_block_ttest.csv', 'wb')
ttest_writer = csv.writer(ttest_file, delimiter=',',
                          quoting=csv.QUOTE_NONE)
ttest_writer.writerow(['gene_expression', 'block', 'p_value'])
for key, gene_block_exp in gene_expression_block.items():
    gene_nonblock_exp = gene_expression_nonblock[key]
    print key
    t_value, p_value = ttest_ind(gene_block_exp, gene_nonblock_exp,
                                 equal_var=False)
    print p_value
    gene_type = key.split('|')[0]
    block_type = key.split('|')[1]
    ttest_writer.writerow([gene_type, block_type, p_value])

block_indi_file.close()
ttest_file.close()


# correlation between methylation difference and gene expression difference

methy_diff = dict()
for tissue_type in tissue_types:
    methy_diff[tissue_type] = []
methy_diff_by_ind = dict()

gene_start_seq = expression_data['start']
gene_end_seq = expression_data['end']
methy_start_seq = methy_raw['start']
methy_end_seq = methy_raw['end']
gene_ind = 0
methy_ind = 0

while gene_ind < gene_start_seq.shape[0]:
    methy_start = methy_start_seq[methy_ind]
    methy_end = methy_end_seq[methy_ind]
    gene_start = gene_start_seq[gene_ind]
    gene_end = gene_end_seq[gene_ind]

    # methy is before gene expression
    if gene_start > methy_end:
        methy_ind += 1
    # methy is after gene expression
    elif gene_end < methy_start:
        gene_ind += 1
    # methy is in between
    elif gene_start <= methy_end <= gene_end or \
        gene_start <= methy_start <= gene_end or methy_start < gene_start < \
            methy_end or methy_start < gene_end < methy_end:
        if gene_ind not in methy_diff_by_ind:
            methy_diff_by_ind[gene_ind] = []
        methy_diff_by_ind[gene_ind].append(methy_ind)
        methy_ind += 1

gene_ind = 0
while gene_ind < gene_start_seq.shape[0]:
    if gene_ind not in methy_diff_by_ind:
        for tissue_type in tissue_types:
            methy_diff[tissue_type].append(0)
    else:
        for tissue_type in tissue_types:
            methy_diff[tissue_type].append(np.mean(methy_raw[tissue_type]
                                                   [methy_diff_by_ind[gene_ind]]))
    gene_ind += 1

print len(methy_diff['breast'])

methy_file = open('../data/gene_methy_correlation.csv', 'wb')
methy_writer = csv.writer(methy_file, delimiter=',',
                          quoting=csv.QUOTE_NONE)
methy_writer.writerow(['gene_expression', 'methylation diff type',
                       'correlation'])

for tissue_type in tissue_types:
    for key, methy_value in methy_diff.items():
        # use the difference of types (normal and tumor) for the same tissue
        expression_type1 = tissue_type + '___normal'
        expression_type2 = tissue_type + '___tumor'
        expression_diff = expression_data[expression_type1] - \
                          expression_data[expression_type2]

        print tissue_type, key
        correlation_coefficient = pearsonr(methy_value,
                                           expression_diff)
        print correlation_coefficient[0]
        methy_writer.writerow([tissue_type, key, correlation_coefficient[
            0]])

methy_file.close()
