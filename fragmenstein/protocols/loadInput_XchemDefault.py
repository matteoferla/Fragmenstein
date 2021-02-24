import os
from collections import OrderedDict

from fragmenstein.protocols.loadInput_base import LoadInput_base
from fragmenstein.utils.compound import Compound
from fragmenstein.utils.io_utils import load_files_as_mols


class LoadInput_XchemDefault(LoadInput_base):

    @staticmethod
    def get_examples_init_params():
        return dict( data_dir= os.path.abspath( os.path.join(__file__, "../../mpro/data/xchem_examples/aligned") ) ,
                     fragIds_pattern=r".*?(x[\w-]+)\.mol$")

    def __init__(self, data_dir, fragIds_pattern, target_id = None):

        super().__init__()
        self.data_dir = data_dir
        self.fragIds_pattern = fragIds_pattern
        self.target_id = target_id



    def prepare_fragments(self):
        fragments = load_files_as_mols( self.data_dir, self.fragIds_pattern )

        def prepareMol(id_mol):
            compound = Compound(id_mol[1])
            compound.molId = id_mol[0]
            return id_mol[0], compound
        fragments = list( map(prepareMol, fragments ) )

        return OrderedDict(fragments)

    def find_templates(self):
       raise NotImplementedError()


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