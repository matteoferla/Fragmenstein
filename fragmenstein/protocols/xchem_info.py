
class Xchem_info():

    fragment_id_pattern = r".*?(x[\w-]+)\.mol$"
    unboundPdb_id_pattern = r".*?(x[\w-]+)_bound\.pdb$"
    unboundPdb_id_pattern = r".*?(x[\w-]+)_apo-desolv\.pdb$"
    @staticmethod
    def default_params_xchem():
        return  dict(

            fragment_id_pattern= Xchem_info.fragment_id_pattern,
            boundPdb_id_pattern= Xchem_info.unboundPdb_id_pattern,
            unboundPdb_id_pattern = Xchem_info.unboundPdb_id_pattern
        )