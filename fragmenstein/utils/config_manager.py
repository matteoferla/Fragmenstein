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

    def parseBool( x ):
        if isinstance(x, str):
            if x == "True":
                return True
            elif x == "False":
                return False
            else:
                raise  ValueError("Error, x (%s) should be either True or False"%x)
        else:
            return bool(x)

    victorVerbose_default =  False if os.environ.get("VICTOR_VERBOSE", None) is None else parseBool(os.environ.get("VICTOR_VERBOSE", None))
    protocolVerbose_default =  False if os.environ.get("PROTOCOL_VERBOSE", None) is None else parseBool(os.environ.get("VICTOR_VERBOSE", None))

    PARSEABLE_PARAMETERS = [
        # parameter_name, str_to_value_funcion, default_value, help
        ("N_CPUS", int, 1, "Number of cpus for parallel computing. Default %(default)s"),
        ("DASK_WORKER_MEMORY", str, "-1", "Memory for each dask worker. E.g 16GB'. Default 90%% of host memory / N_CPUS"),
        ("VICTOR_VERBOSE", parseBool, victorVerbose_default, "Enable verbosity in Victor. Default: %(default)s"),
        ("PROTOCOL_VERBOSE", parseBool, protocolVerbose_default, "Enable verbosity in protocol. Default: %(default)s"),
        ("TMP_DIR", str, "/tmp", "Temporal directory for computations. Default: %(default)s"),
        ("RANDOM_SEED_PERMUT", int, 121, "Random seed for permutations. Default: %(default)s"),
    ]

    NONPARSEABLE_PARAMETERS = [
        ("EXTERNAL_TOOLS_CONFIG_FILE", str, os.path.abspath(os.path.join(__file__, "../../external/external_config.json")),  "A json file for external tools configuration"),
        ("COMBINE_PERMUTATIONS_MAX_NUM_ELEMENTS", int, 2, "The max number of fragments that would be considered when doing permutaitons of fragments combinations. Default: %(default)s"),
        ("COMBINE_PERMUTATIONS_MAX_NUM_ELEMENTS", int, 2,
         "The max number of fragments that would be considered when doing permutaitons of fragments combinations. Default: %(default)s"),
        ("DISABLE_DASHBOARD", parseBool, False, "Set it to True if you want to disable dask dashboard. Default: %(default)s"),
        #("IGOR_TIMEOUT", int, 200, "Seconds till rosetta minimize timeout. Default: %(default)s") #Mot working because of dask nanny
        ("DE_LINKER_N_ATOMPAIRS_ATTEMPS", int, 3,  "Number of minimun distance pairs that will be intended to be linked"),
    ]
    def __init__(self):
        self_class = type(self)
        self.parseable_params = list([elem[0] for elem in self_class.PARSEABLE_PARAMETERS])
        names, funcs, defaults, helps = zip(* (self_class.PARSEABLE_PARAMETERS+self_class.NONPARSEABLE_PARAMETERS) )
        self.helps = helps
        self.params_dict = dict(zip( names, funcs) )
        self.default_params = dict(zip( names, defaults) )

    def _isBoolean(self, name):
        return isinstance(self.default_params[name], bool)

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

    def yieldParseableParams(self):
        for name, help in zip(ConfigManager.default_params.keys(), ConfigManager.helps):
            if name in self.parseable_params:
                yield name, help


if __name__ == "__main__":

    print( ConfigManager.TMP_DIR)
    print( ConfigManager.N_CPUS)
    for key, val in ConfigManager.default_params.items():
        print(key, val)