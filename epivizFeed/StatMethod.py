class StatMethod(object):
    def __init__(self, measurements=None):
        self.measurements = measurements

    def compute(self, chromosome, start, end):
        print("base class")

    def get_measurements_self(self, mes_type):
        measurement_types = []
        for measurement in self.measurements:
            self.measurement_type_switch(measurement_types, measurement, mes_type)

        return measurement_types

    def measurement_type_switch(self, measurement_types, measurement, mes_type):
        data_obj = {
                        "id": measurement["id"],
                        "name": measurement["name"],
                        "datasourceId": measurement["datasourceId"]
                    }
        if mes_type == "gene":
            if measurement["defaultChartType"] == "scatterplot":
                measurement_types.append(data_obj)

        elif mes_type == "block":
            if measurement["defaultChartType"] == "block":
                measurement_types.append(data_obj)

        elif mes_type == "methy":
            if measurement["defaultChartType"] == "line":
                if measurement["datasourceId"] == "timp2014_probelevel_beta":
                    measurement_types.append(data_obj)

        elif mes_type == "methy_diff":
            if measurement["defaultChartType"] == "line":
                if measurement["datasourceId"] != "timp2014_probelevel_beta":
                    measurement_types.append(data_obj)
