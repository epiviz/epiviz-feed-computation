from urllib.request import urlopen
from requests import get_methy_data

class methy:
	def __init__(self, start, end, chromesome, measurement):
		if measurement["defaultChartType"] == "line":
						raise ValueError('nonnonononno')
		else if measurement["datasourceId"] == "timp2014_probelevel_beta":
			self.methy = get_methy_data(start,end,chromosome, measurement)
		else 
			self.diff_methy = get_methy_data(start,end,chromosome, measurement)
	def calculate_methy_diff():
        methy_corr_res = []
        # loop through every possible combinations of methylation
        for data_source_one, data_source_two in itertools.combinations(
                methylation_diff_types, 2):
            type1 = data_source_one["id"]
            type2 = data_source_two["id"]
            correlation_coefficient = pearsonr(methy_raw_diff[type1], methy_raw_diff[
                type2])
            data_range = {
                'attr-one': [min(methy_raw_diff[type1]), max(methy_raw_diff[type1])],
                'attr-two': [min(methy_raw_diff[type2]), max(methy_raw_diff[type2])]
            }
            corr_obj = build_obj('correlation', 'methylation diff',
                                 'methylation diff', True, data_source_one,
                                 data_source_two,
                                 correlation_coefficient[0],
                                 correlation_coefficient[1],
                                 ranges=data_range)
            methy_corr_res.append(corr_obj)
        methy_corr_res = sorted(methy_corr_res, key=lambda x: x['value'],
                                reverse=True)
        yield methy_corr_res
    
        def methy_res():
    	methy_corr_res = []
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