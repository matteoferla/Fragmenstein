import os
from rdkit import Chem
from fragmenstein import Victor

from fragmenstein.external import ExternalToolImporter
import logging

[pyrosetta]= ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])

project_dir = os.path.expanduser("~/tmp/frament_network")
smiles_fname= os.path.join(project_dir, "input/x0107_x1382_s2.smi")
pdb_filename = os.path.join(project_dir, "input/Mpro-x0107_0A_apo-desolv.pdb")
fragments_fnames = [ os.path.join(project_dir, "input", basename) for basename in ["Mpro-x0107_0A.mol", "Mpro-x1382_0A.mol"] ]
out_dir = os.path.join(project_dir, "output")



hits = [ Chem.MolFromMolFile(frag) for frag in fragments_fnames ] #A list of fragments

Victor.enable_stdout(level=logging.DEBUG)
#Victor.error_to_catch = NotImplementedError #Uncomment if you want fragmenstein to break on error.


Victor.work_path=  out_dir #where results would be saved

v= Victor(       hits=hits,
                 pdb_filename= pdb_filename, #file name of apo protein
                 ligand_resn="LIG",
                 covalent_resn= 'CYS', covalent_resi= '145A'
)


smiles = []
with open( smiles_fname) as f:
    for i, line in enumerate(f):
        # if i < 72: continue
        smi = line.strip()
        smiles.append( smi )


v.place(smiles=smi,
    long_name= '%d'%i)

v.make_pse()
print( v.summarise() )


#TODO: use the protocols instead

with open( smiles_fname) as f:
    for i, line in enumerate(f):
        if i < 72: continue
        smi = line.strip()
        print(smi)
        v.place(smiles=smi,
            long_name= '%d'%i)

        v.make_pse()
        print( v.summarise() )


'''
python -m examples.fragNet_followup_fragmenstein
'''