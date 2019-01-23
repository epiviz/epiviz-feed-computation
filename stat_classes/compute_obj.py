from stat_classes.ttest_block import Ttest_Block
from stat_classes.ttest_gene import Ttest_Gene
from stat_classes.block_overlap_percent import Block_overlap


class Compute_obj:

    def init(self, computation, measurements):
        self.computation_obj = self.create(computation)
        self.measurements = measurements

    def create(self, computation):

        if(computation == "ttest_block_expression"):
            ret_val = Ttest_Block(self.measurements)
        elif(computation == "block_overlap_percent"):
            ret_val = Block_overlap(self.measurements)
        elif(computation == "expression_methydiff_correlation"):
            ret_val = expression_methydiff_correlation(self.measurements)
        elif(computation == "expression_methy_correlation"):
            ret_val = ttest_obj
        elif(computation == "ttest_expression_per_gene"):
            ret_val = Ttest_Gene(self.measurements)
        elif(computation == "methy_diff_correlation"):
            ret_val = ttest_obj
        elif(computation == "expression_correlation"):
            ret_val = ttest_obj

        return ret_val

    def compute(self, chromosome, start_seq, end_seq):
        if(self.computation == "ttest_block_expression"):
            ret_val = self.computation_obj.compute()
        elif(self.computation == "block_overlap_percent"):
            ret_val = self.computation_obj.compute(start_seq, end_seq)
        elif(self.computation == "expression_methydiff_correlation"):
            ret_val = self.computation_obj.compute()
        elif(self.computation == "expression_methy_correlation"):
            ret_val = self.computation_obj.compute()
        elif(self.computation == "ttest_expression_per_gene"):
            ret_val = self.computation_obj.compute(chromosome, start_seq, end_seq)
        elif(self.computation == "methy_diff_correlation"):
            ret_val = self.computation_obj.compute()
        elif(self.computation == "expression_correlation"):
            ret_val = self.computation_obj.compute()

        return ret_val
