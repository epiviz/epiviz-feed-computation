
def computation_request(start_seq, end_seq, chromosome, gene_name, measurements=None):
    # extract data from measurements
    gene_types = []
    block_types = []

    methylation_types = []
    methylation_diff_types = []
    # categorize measurements into different types
    for measurement in measurements:
        data_obj = {
            "id": measurement["id"],
            "name": measurement["name"],
            "datasourceId": measurement["datasourceId"]
        }
        if measurement["defaultChartType"] == "scatterplot":
            gene_types.append(data_obj)
        elif measurement["defaultChartType"] == "block":
            block_types.append(data_obj)
        elif measurement["defaultChartType"] == "line":
            if measurement["datasourceId"] == "timp2014_probelevel_beta":
                methylation_types.append(data_obj)
            else:
                methylation_diff_types.append(data_obj)

    block_data = None
    methy_raw = None
    methy_raw_diff = None
    expression_data = None
    has_block = len(block_types) > 0
    has_methy = len(methylation_types) > 0
    has_methy_diff = len(methylation_diff_types) > 0
    has_gene = len(gene_types) > 0

    expression_data = get_gene_data(start_seq, end_seq, chromosome,
                                    gene_types)

    per_gene_ttest = ttest_expression_per_gene(gene_types, expression_data,
                                               chromosome, start_seq, end_seq)

    yield per_gene_ttest

    if has_block:
        block_data = get_block_data(
            start_seq, end_seq, chromosome, block_types)

        # block overlap percentage
        block_overlap = block_overlap_percent(block_types, block_data,
                                              start_seq, end_seq)
        yield block_overlap

    if has_methy_diff:
        methy_raw_diff = get_methy_data(start_seq, end_seq, chromosome,
                                        methylation_diff_types)

        methy_diff_corr_res = methy_diff_correlation(
            methy_raw_diff, methylation_diff_types)

        yield methy_diff_corr_res

    if has_methy:
        methy_raw = get_methy_data(start_seq, end_seq, chromosome,
                                   methylation_types)
        methy_corr_res = methy_correlation(methy_raw, methylation_diff_types)

        # loop through normal/tumor of each tissue type
        for data_source_one, data_source_two in itertools.combinations(
                methylation_diff_types, 2):
            type1 = data_source_one["id"]
            type2 = data_source_two["id"]
            if type1.split("_")[0] != type2.split("_")[0]:
                continue

            correlation_coefficient = pearsonr(methy_raw[type1], methy_raw[
                type2])

            print (type1, type2)
            data_range = {
                'attr-one': [min(methy_raw[type1]), max(methy_raw[type1])],
                'attr-two': [min(methy_raw[type2]), max(methy_raw[type2])]
            }
            corr_obj = build_obj('correlation', 'methylation',
                                 'methylation', True, data_source_one,
                                 data_source_two,
                                 correlation_coefficient[0],
                                 correlation_coefficient[1],
                                 ranges=data_range)
            methy_corr_res.append(corr_obj)
        methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
                                reverse=True)
        yield methy_corr_res

    if has_gene:
        # expression_data = get_gene_data(start_seq, end_seq, chromosome,
        #                                 gene_types)

        corr_list = []
        # pvalue_list = []
        for data_source_one, data_source_two in itertools.combinations(
                gene_types, 2):
            exp1 = data_source_one['id']
            exp2 = data_source_two['id']

            if exp1 not in expression_data.columns or exp2 not in expression_data.columns:
                continue

            col_one = expression_data[exp1]
            col_two = expression_data[exp2]

            correlation_coefficient = pearsonr(col_one, col_two)
            corr_obj = build_obj('correlation', 'expression', 'expression',
                                 True, data_source_one,
                                 data_source_two, correlation_coefficient[0],
                                 correlation_coefficient[1])
            corr_list.append(corr_obj)

            t_value, p_value = ttest_ind(col_one, col_two,
                                         equal_var=False)
            # ttest_obj = build_obj('t-test', 'expression', 'expression', True,
            #                       data_source_one, data_source_two, t_value,
            #                       p_value)
            # pvalue_list.append(ttest_obj)

        # pvalue_list = sorted(pvalue_list, key=lambda x: x['value'],
        #                      reverse=True)
        # yield pvalue_list

        corr_list = sorted(corr_list, key=lambda x: x['value'],
                           reverse=True)
        yield corr_list

    if has_gene and has_block:
        # gene expression and block independency test
        ttest_block_exp = ttest_block_expression(expression_data, block_data,
                                                 gene_types, block_types)
        yield ttest_block_exp

    if has_gene and has_methy:
        # correlation between methylation and gene expression
        # with the same tissue type
        corr_methy_gene = expression_methy_correlation(expression_data, gene_types, methylation_types,
                                                       methy_raw)

        yield corr_methy_gene

    if has_gene and has_methy_diff:
        # correlation between methylation difference and gene expression
        # difference
        corr_methy_gene = expression_methydiff_correlation(expression_data, gene_types, methylation_diff_types,
                                                           methy_raw_diff)

        yield corr_methy_gene
