
def format_exp_methy_output(attr1, attr2, type1, type2):
    data = []
    for exp, methy in zip(attr1, attr2):
        point = dict()
        point[type1] = exp
        point[type2] = methy
        data.append(point)

    return data
