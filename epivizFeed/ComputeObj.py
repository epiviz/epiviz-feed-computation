from stat_classes.TtestBlock import TtestBlock
from stat_classes.TtestGene import TtestGene
from stat_classes.OverlapBlockPercent import OverlapBlock
from stat_classes.CorrelationExpMethy import CorrelationExpMethy
from stat_classes.CorrelationExp import CorrelationExp
from stat_classes.CorrelationMethy import CorrelationMethy


class ComputeObj:

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
        elif(computation == "methy_diff_correlation"):
            ret_val = CorrelationMethy(self.measurements, "methy_diff")
        elif(computation == "methy_correlation"):
            ret_val = CorrelationMethy(self.measurements, "methy")
        elif(computation == "expression_correlation"):
            ret_val = CorrelationExp(self.measurements)

        return ret_val

    def compute(self, chromosome, start, end, additional=None):

        return self.computation_obj.compute(chromosome, start, end)
