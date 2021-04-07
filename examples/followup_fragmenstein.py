import os
from rdkit import Chem
from fragmenstein import MOLS_EXAMPLES_DIR, Victor

from fragmenstein.external import ExternalToolImporter
import logging

[pyrosetta]= ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])

smiles_fname= os.path.join(MOLS_EXAMPLES_DIR, "followup_x20706.smiles")
pdb_filename = os.path.join(MOLS_EXAMPLES_DIR, "F709-PHIPA-x20706_unbound.pdb")  #Cleaned protein with no water, no ligand.
fragments_fnames = [ os.path.join(MOLS_EXAMPLES_DIR, "F709-PHIPA-x2070.mol") ]
out_dir = "ouptut" # os.path.expanduser("~/tmp/output_followup")

with open( smiles_fname) as f:
  smi = f.read().strip()
  print(smi)

hits = [ Chem.MolFromMolFile(frag) for frag in fragments_fnames ] #A list of fragments

Victor.enable_stdout(level=logging.DEBUG)
#Victor.error_to_catch = NotImplementedError #Uncomment if you want fragmenstein to break on error.


Victor.work_path=  out_dir #where results would be saved

v= Victor(       hits=hits,
                 pdb_filename= pdb_filename, #file name of apo protein
                 ligand_resn="LIG",
                 # #Next line is nonsense for non covalent residues but required (at the moment). Just pick a random CYS
                 # covalent_resn= 'CYS',
                  covalent_resi= '1382A'
)

v.place(smiles=smi,
        long_name= 'x20706_out', #just to name output files
)

v.make_pse()
print( v.summarise() )


'''
python -m examples.followup_fragmenstein
'''