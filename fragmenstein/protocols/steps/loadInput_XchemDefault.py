import os
import re
from collections import OrderedDict

from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.loadInput_base import LoadInput_base
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.utils.io_utils import load_files_as_mols, apply_func_to_files


class LoadInput_XchemDefault(LoadInput_base, Xchem_info):

    @staticmethod
    def get_examples_init_params():
        params = dict( data_dir= os.path.abspath( os.path.join(__file__, "../../mpro/data/xchem_examples/aligned") ) )
        params.update( Xchem_info.default_params_xchem())
        return params

    def __init__(self, data_dir, fragment_id_pattern, unboundPdb_id_pattern, target_id = None, *args, **kwargs):

        super().__init__()
        self.data_dir = data_dir
        self.fragment_id_pattern = fragment_id_pattern
        self.unboundPdb_id_pattern = unboundPdb_id_pattern
        self.target_id = target_id



    def prepare_fragments(self):
        fragments = load_files_as_mols( self.data_dir, self.fragment_id_pattern )

        def prepareMol(id_mol):
            compound = Compound(id_mol[1])
            compound.molId = id_mol[0]
            return id_mol[0], compound
        fragments = list( map(prepareMol, fragments ) )

        return OrderedDict(fragments)

    def find_templates(self, filter_ids=None):
        templates = apply_func_to_files(self.data_dir, self.unboundPdb_id_pattern,
                                        lambda x: (re.match(self.unboundPdb_id_pattern,x).group(1), x) )
        if filter_ids:
            templates = list(filter(lambda id_fname: id_fname[0] in filter_ids, templates))

        return templates


def test_load():


    lxchem = LoadInput_XchemDefault( ** LoadInput_XchemDefault.get_examples_init_params() )
    print(next(lxchem.fragments))
    print( lxchem.fragments_dict)
    print(lxchem.fragments_dict["x2646_0A"].molId)
    pass

if __name__ =="__main__":
    test_load()

    '''

python -m fragmenstein.protocols.loadInput_XchemDefault

    '''