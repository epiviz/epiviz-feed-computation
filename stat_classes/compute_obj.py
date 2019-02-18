from stat_classes.ttest_block import TtestBlock
from stat_classes.ttest_gene import TtestGene
from stat_classes.overlap_block_percent import OverlapBlock
from stat_classes.correlation_exp_methy import CorrelationExpMethy


class Compute_obj:

    def __init__(self, computation, measurements):
        self.measurements = measurements
        self.computation_obj = self.create(computation)

    def create(self, computation):

        if(computation == "ttest_block_expression"):
            ret_val = TtestBlock(self.measurements)
        elif(computation == "block_overlap_percent"):
            ret_val = OverlapBlock(self.measurements)
        elif(computation == "expression_methydiff_correlation"):
            ret_val = CorrelationExpMethy(self.measurements, "methy_diff")
        elif(computation == "expression_methy_correlation"):
            ret_val = CorrelationExpMethy(self.measurements, "methy")
        elif(computation == "ttest_expression_per_gene"):
            ret_val = TtestGene(self.measurements)
        # elif(computation == "methy_diff_correlation"):
        #     ret_val = ttest_obj
        # elif(computation == "expression_correlation"):
        #     ret_val = ttest_obj

        return ret_val

    def compute(self, chromosome, start, end):

        return self.computation_obj.compute(chromosome, start, end)