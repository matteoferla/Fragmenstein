import sys, os
import pyrosetta
from rdkit import Chem

pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

sys.path.append("/home/ruben/oxford/tools/Fragmenstein")

from fragmenstein import Monster, Victor, Igor, Rectifier


smiles_fname="./followup_x20706.smiles"
pdb_filename ="./F709-PHIPA-x20706_unbound.pdb" #Cleaned protein with no water, no ligand.

with open( smiles_fname) as f:
  smiles = f.read().strip()
  print( smiles )

hits = [ Chem.MolFromMolFile("./F709-PHIPA-x2070.mol") ] #A list of fragments

import logging
Victor.enable_stdout(level=logging.DEBUG)
#Victor.error_to_catch = NotImplementedError #Uncomment if you want fragmenstein to break on error.

Victor.output_dir= "./" #where results would be saved
v= Victor(smiles=smiles, hits=hits,
                 pdb_filename= pdb_filename, #file name of apo protein
                 long_name= 'x20706_out', #just to name output files
                 #Next line is nonsense for non covalent residues but required (at the moment). Just pick a random CYS
                 covalent_resn= 'CYS', covalent_resi= '1382A', 
                 )
v.make_pse()
print( v.summarise() )

