from stat_classes.compute_obj import Compute_obj
from run import test_measurements

measurements = test_measurements()
computations = ["ttest_block_expression", "block_overlap_percent",
                "expression_methydiff_correlation", "expression_methy_correlation",
                "ttest_expression_per_gene", "methy_diff_correlation", "methy_correlation",
                "expression_correlation"]
chromosome, start, end = "chr1", 10550488, 11554489
for comp in computations:
    obj = Compute_obj(comp, measurements)
    result = obj.compute(chromosome, start, end)
    print(result)
