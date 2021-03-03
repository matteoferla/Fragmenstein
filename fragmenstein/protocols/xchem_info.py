
class Xchem_info():

    fragment_id_pattern = r".*?(x[\w-]+)\.mol$"
    boundPdb_id_pattern = r".*?(x[\w-]+)_bound\.pdb$"
    unboundPdb_id_pattern = r".*?(x[\w-]+)_apo-desolv\.pdb$"
    predicted_boundPdb_id_pattern = r".*?(x[\w\-_]+)\.holo_minimised\.pdb$"
    fragment_no_chain_pattern = r"^.*(x\d+).*$"
    @staticmethod
    def default_params_xchem():
        return  dict(

            fragment_id_pattern= Xchem_info.fragment_id_pattern,
            boundPdb_id_pattern= Xchem_info.boundPdb_id_pattern,
            unboundPdb_id_pattern = Xchem_info.unboundPdb_id_pattern,
            predicted_boundPdb_id_pattern = Xchem_info.predicted_boundPdb_id_pattern
        )