import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

from fragmenstein.scoring.scoring_config import MPRO_HITS_DIR

BAD_RESIDUES_LIST = ["DMS"]
LIGAND_TAG = "LIG"

def load_hits(hits_folder= MPRO_HITS_DIR):
    hits_dict = {}
    for fname in os.listdir(hits_folder):
        hit_code = fname.split("-")[1].split(".")[0]
        hits_dict[hit_code] = Chem.MolFromMolFile( os.path.join(hits_folder, fname) )
    return hits_dict


fname_pdbBound = os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x10355_0A/Mpro-x10355_0A_bound.pdb")
pdb_mol = Chem.MolFromPDBFile( fname_pdbBound )


bad_residues =  Chem.rdmolops.SplitMolByPDBResidues(pdb_mol, whiteList= BAD_RESIDUES_LIST, negateList=False )
for br in bad_residues.values():
    pdb_mol = Chem.rdmolops.DeleteSubstructs(pdb_mol, br )

ligand = Chem.rdmolops.SplitMolByPDBResidues(pdb_mol, whiteList= [LIGAND_TAG] )[LIGAND_TAG]
protein = Chem.rdmolops.DeleteSubstructs(pdb_mol, ligand )

protein_conf = protein.GetConformer()

print("Protrusion ligand-protein: %2.4f"%Chem.rdShapeHelpers.ShapeProtrudeDist(ligand, protein) )
print("Protrusion ligand-ligand: %2.4f"%Chem.rdShapeHelpers.ShapeProtrudeDist(ligand, ligand))
print("Protrusion protein-protein: %2.4f"%Chem.rdShapeHelpers.ShapeProtrudeDist(protein, protein))


hits_dict = load_hits()

protrusions_list = []
for hit_id, hit_mol in hits_dict.items():
    protrusion = Chem.rdShapeHelpers.ShapeProtrudeDist(ligand, hit_mol) / Descriptors.ExactMolWt(hit_mol)
    protrusions_list.append( (hit_id, protrusion) )

print( sorted(protrusions_list, key = lambda x: x[-1], reverse = True ) )

'''
python -m fragmenstein.scoring.volumeTrials
'''