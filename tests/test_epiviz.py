import sys
sys.path.insert(0,'C:/Users/Kyle/Documents/Python Scripts/Envs/research_epiviz/efeedcomputation/src/epivizfeedcompute')

from stats import BaseStats, TtestBlock, TtestExp, Correlation, CorrelationGeneSignal, OverlapBlock
# __author__ = "Jayaram Kancherla"
# __copyright__ = "Jayaram Kancherla"
# __license__ = "mit"
overlap_measurements = [
    {
        "id": "breast___normal",
        "name": "Expression breast_normal",
        "type": "feature",
        "datasourceId": "gene_expression_barcode_subtype",
        "datasourceGroup": "gene_expression_barcode_subtype",
        "dataprovider": "umd",
        "formula": None,
        "datatype": "expression",
        "defaultChartType": "scatterplot",
        "annotation": None,
        "metadata": ["probe"]
    },
    {
        "id": "breast___tumor",
        "name": "Expression breast_tumor",
        "type": "feature",
        "datasourceId": "gene_expression_barcode_subtype",
        "datasourceGroup": "gene_expression_barcode_subtype",
        "dataprovider": "umd",
        "formula": None,
        "datatype": "expression",
        "defaultChartType": "scatterplot",
        "annotation": None,
        "metadata": ["probe"]
    }
]

chrom = "chr11"
start = 10000000
end = 1235623
def test_ttest():
    # create instance of the class
    test = Correlation.Correlation(overlap_measurements, 0.05)
    result = test.compute(chrom, start, end, {"datatype": "expression", "annotation":None})
    print(result)
    assert result

def test_overlap():
    # create instance of the class
    test = OverlapBlock.OverlapBlock(overlap_measurements, 0.05)
    result = test.compute(chrom, start, end, {"datatype": "expression", "annotation":None})
    print(result)
    assert result
    
test_ttest()