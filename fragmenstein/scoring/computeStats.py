import os
import tempfile

import pyrosetta
import sys
from rdkit import Chem
from rdkit_to_params import Params

from fragmenstein import Igor
from fragmenstein import Victor
from fragmenstein.external import ExternalToolImporter
from fragmenstein.mpro import MProVictor
from fragmenstein.scoring.scoring_config import MPRO_RAW_DATA_DIR

compound_xid = "x1336"
fragments = [ "x1478", "x1493"]

all_compounds_aligned_path = os.listdir( os.path.join(MPRO_RAW_DATA_DIR, "aligned") )

compoundName = [ elem for elem in all_compounds_aligned_path if elem.startswith("Mpro-"+compound_xid)][0]

fullPrefix = os.path.join(MPRO_RAW_DATA_DIR, "aligned", compoundName)

molfile = os.path.join(fullPrefix , compoundName+".mol")
mol = Chem.MolFromMolFile(molfile)
print( Chem.MolToSmiles(mol ) )

bound_pdb_fname= os.path.join(fullPrefix, compoundName+"_bound.pdb")
apo_pdb_fname= os.path.join(fullPrefix, compoundName+"_apo-desolv.pdb")

mPro_pdb = os.path.abspath(os.path.join(MProVictor.get_mpro_path(), "template.pdb"))


smiles= Chem.MolToSmiles(mol)

# pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
#
# MProVictor.quick_renanimation = True
# victor = MProVictor.from_hit_codes(smiles= Chem.MolToSmiles(mol),
#                                    hit_codes=fragments,
#                                    long_name=compoundName)
# print( victor.summarise() )


[sascorer]=  ExternalToolImporter.import_tool("DeLinker", ["sascorer"])


sa = sascorer.calculateScore( mol )
print(sa)

'''
    LE = (ΔG)/N
'''