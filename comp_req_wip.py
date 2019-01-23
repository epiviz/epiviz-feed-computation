from stat_classes.compute_obj import Compute_obj


def computational_request(config):
    measurements, computations = init(config)
    chromosome, start, end = regions_of_interest()
    for comp in computations:
        obj = Compute_obj(comp, measurements)
        result = obj.compute(chromosome, start, end)
        yield result
