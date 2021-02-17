import os

#TOOD: Implement a cmd parser over argparse parser.parse_known_args()

def instancer(cls):
    '''
    Will create an instance of a given class
    :param cls: class
    :return: instance of a class
    '''
    return cls()

@instancer #ConfigManager will be an instance of the ConfigManager() class
class ConfigManager():

    #ADD any parameter you want to PARAMETERS. #TODO: Load default values from some config file
    PARAMETERS = [
        # parameter_name, str_to_value_funcion, default_value, help
        ("N_CPUS", int, 1, "Number of cpus for parallel computing. Default %(default)s"),
        ("DASK_WORKER_MEMORY", str, None, "Memory for each dask worker. E.g '16GB"),
        ("TMP_DIR", str, "/tmp", "Temporal directory for computations %(default)s"),
    ]

    def __init__(self):
        names, funcs, defaults, helps = zip(* type(self).PARAMETERS)
        self.helps = helps
        self.params_dict = dict(zip( names, funcs) )
        self.default_params = dict(zip( names, defaults) )

    def get(self, name):
        assert  name in self.params_dict, "Error, {name} is not a configurable parameter".format(name=name)
        property = os.environ.get(name, None)
        if property:
            property = self.params_dict[name](property)
        else:
            property = self.default_params[name]
        return property

    def __getattr__(self, name):
        return self.get(name)

if __name__ == "__main__":

    print( ConfigManager.TMP_DIR)
    print( ConfigManager.N_CPUS)