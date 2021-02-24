
class Xchem_info():

    fragment_id_pattern = r".*?(x[\w-]+)\.mol$"
    boundPdb_id_pattern = r".*?(x[\w-]+)_bound\.pdb$"

    @staticmethod
    def default_params_xchem():
        return  dict(

            fragment_id_pattern= Xchem_info.fragment_id_pattern,
            boundPdb_id_pattern= Xchem_info.boundPdb_id_pattern
        )