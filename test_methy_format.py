# test file for getting average methylation data according to gene expression
# regions

from methy_format_gene import methy_format_gene

test_measurements = [{
    'name': 'expression colon',
    'id': 'colon___normal'
}, {
    'name': 'expression breast',
    'id': 'breast___tumor'
}, {
    'name': 'methylation colon',
    'id': 'colon'
}, {
    'name': 'methylation lung',
    'id': 'lung'
}]

results = methy_format_gene(3947953, 7164991, 'chr11', test_measurements)

print results
