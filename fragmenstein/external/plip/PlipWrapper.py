import os
from collections import defaultdict

import sys

from fragmenstein.external import ExternalToolImporter

[plip_preparation] =  ExternalToolImporter.import_tool("plip",
                                       ["plip.structure.preparation"])

class PlipWrapper():

    AVALIABLE_INTERACTIONS = ['saltbridge_lneg', 'saltbridge_pneg','hbonds_ldon', 'hbonds_pdon', 'pistacking','pication_paro', 'pication_laro',
                              'hydrophobic_contacts', 'halogen_bonds', 'water_bridges', 'metal_complexes']
    AVAILABLE_UNPAIRED = ['unpaired_hba', 'unpaired_hbd', 'unpaired_hal']

    def __init__(self, ligand_resname: str="LIG"):
        '''
        :param ligand_resname: str. The residue name of the ligand.
        '''
        self.ligand_resname = ligand_resname

    def compute_interactions_boundPDB(self, fname):
        '''
        :param fname: str. a filename
        '''
        my_mol = plip_preparation.PDBComplex()
        my_mol.load_pdb( fname)
        my_mol.analyze()

        interactions = [ my_mol.interaction_sets[keyId] for keyId in my_mol.interaction_sets if keyId.startswith(self.ligand_resname) ]
        interactions_dict = defaultdict(list)
        for inters_list in interactions:
            for inter_type in PlipWrapper.AVALIABLE_INTERACTIONS:
                interaction_instances = getattr(inters_list, inter_type)
                for inter in interaction_instances:
                    binding_id = "%s_%s_%s"%(inter.resnr, inter.restype,  inter.reschain)
                    interactions_dict[ inter_type ].append( binding_id )
        return interactions_dict

if __name__ == "__main__":

    fname = os.path.expanduser( sys.argv[1])

    pw = PlipWrapper()
    results = pw.compute_interactions_boundPDB(fname)
    print( list(results.items()))

    '''
python -m fragmenstein.external.plip.PlipWrapper ~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x3124_0A/Mpro-x3124_0A_bound.pdb

    '''