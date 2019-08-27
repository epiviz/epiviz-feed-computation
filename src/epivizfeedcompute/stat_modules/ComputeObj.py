from .TtestBlock import TtestBlock
from .TtestGene import TtestGene
from .OverlapBlockPercent import OverlapBlock
from .CorrelationExpMethy import CorrelationExpMethy
from .CorrelationExp import CorrelationExp
from .CorrelationMethy import CorrelationMethy

class ComputeObj:

    def __init__(self, computation, measurements, pval_threshold):
        self.measurements = measurements
        self.computation_obj = self.create(computation, pval_threshold)

    def create(self, computation, pval_threshold):

        if(computation == "ttest_block_expression"):
            ret_val = TtestBlock(self.measurements, pval_threshold)
        elif(computation == "block_overlap_percent"):
            ret_val = OverlapBlock(self.measurements, pval_threshold)
        elif(computation == "expression_methydiff_correlation"):
            ret_val = CorrelationExpMethy(self.measurements, "methy_diff", pval_threshold)
        elif(computation == "expression_methy_correlation"):
            ret_val = CorrelationExpMethy(self.measurements, "methy", pval_threshold)
        elif(computation == "ttest_expression_per_gene"):
            ret_val = TtestGene(self.measurements, pval_threshold)
        elif(computation == "methy_diff_correlation"):
            ret_val = CorrelationMethy(self.measurements, "methy_diff", pval_threshold)
        elif(computation == "methy_correlation"):
            ret_val = CorrelationMethy(self.measurements, "methy", pval_threshold)
        elif(computation == "expression_correlation"):
            ret_val = CorrelationExp(self.measurements, pval_threshold)

        return ret_val

    def compute(self, chromosome, start, end, additional=None):

        return self.computation_obj.compute(chromosome, start, end)
