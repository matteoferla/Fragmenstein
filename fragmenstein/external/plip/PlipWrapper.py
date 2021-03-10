import os
import tempfile
from collections import defaultdict

import sys
from typing import Dict, List

from fragmenstein.external import ExternalToolImporter



class PlipWrapper():
    '''
    A wrapper class to use plip in a transparent manner. PLIP is a python package that analyzes molecular interactions in pdb files
    '''
    AVALIABLE_INTERACTIONS = ['saltbridge_lneg', 'saltbridge_pneg','hbonds_ldon', 'hbonds_pdon', 'pistacking','pication_paro', 'pication_laro',
                              'hydrophobic_contacts', 'halogen_bonds', 'water_bridges', 'metal_complexes']
    AVAILABLE_UNPAIRED = ['unpaired_hba', 'unpaired_hbd', 'unpaired_hal']

    #importing plip as a class attribute to avoid loading it when importing
    [plip_preparation] = ExternalToolImporter.import_tool("plip",
                                                          ["plip.structure.preparation"])
    def __init__(self, ligand_resname: str="LIG"):
        '''
        :param ligand_resname: str. The residue name of the ligand in bound pdbs
        '''
        self.ligand_resname = ligand_resname

    def compute_interactions_boundPDB(self, fname:str) -> Dict[str, List[str]]:
        '''
        computes the interactions between the protein and the ligand
        :param fname: str. a pdb filename
        :return a dictionary that containes, for each type interaction, the protein residue ids that interact with the ligand
                e.g. [('hbonds_pdon', ['166_A_GLU', '142_A_ASN', '143_A_GLY']), ('hydrophobic_contacts', ['166_A_GLU', '187_A_ASP', '189_A_GLN'])]
        '''


        cur_wd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            my_mol = PlipWrapper.plip_preparation.PDBComplex()
            my_mol.load_pdb( fname)
            my_mol.analyze()
            os.chdir(cur_wd)
        interactions = [ my_mol.interaction_sets[keyId] for keyId in my_mol.interaction_sets if keyId.startswith(self.ligand_resname) ] #TODO. Include distance and angle
        interactions_dict = defaultdict(list)
        for inters_list in interactions:
            for inter_type in PlipWrapper.AVALIABLE_INTERACTIONS:
                interaction_instances = getattr(inters_list, inter_type)
                for inter in interaction_instances:
                    binding_id = "%s_%s_%s"%(inter.resnr, inter.reschain, inter.restype)
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