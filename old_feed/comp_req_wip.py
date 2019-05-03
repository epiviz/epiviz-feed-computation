import pandas as pd
from stat_classes.ComputeObj import ComputeObj
from run import test_measurements

measurements = test_measurements()
computations = ["ttest_block_expression"]

# computations = ["ttest_block_expression", "block_overlap_percent",
#                 "expression_methydiff_correlation", "expression_methy_correlation",
#                 "ttest_expression_per_gene", "methy_diff_correlation", "methy_correlation",
#                 "expression_correlation"]
chromosome, start, end = "chr1", 10550488, 11554489
results = pd.DataFrame()
for comp in computations:
    obj = ComputeObj(comp, measurements)
    results = results.append(obj.compute(chromosome, start, end))
print(results)
