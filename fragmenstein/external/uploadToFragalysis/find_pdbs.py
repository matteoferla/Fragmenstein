import os

import shutil
from rdkit import Chem



'''
ori_name = "x0034-0B-0b0-x0176-0B-0b0-x0176-0B-0b1-x0438-0B-0b2-x0438-0B-0b4-x0438-0B-0b5"
x0034-0B-x0212-0B.holo_unminimised.pdb


'''


wdir_root = "/data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/enumeration/"
wdirs = [ os.path.join(wdir_root, name) for name in ["Site_1_brics", "Site_1_perm"]]

sdf_fname = "/data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/enumeration/final/1B_nsp13_nucletotide_fragmestein_filtered.sdf"
out_dir = "/data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/enumeration/Site_1_fragmestein_pdbs"

suppl = Chem.SDMolSupplier(sdf_fname)

n_found = 0
founds ={}

i = 0
for i, mol in enumerate(suppl):
    if i==0: continue
    print(mol.GetProp("_Name"))
    ori_name = mol.GetProp("original_name")
    for wdir in wdirs:
        pdb_fname = os.path.join(wdir, "merges", ori_name, "%s.holo_minimised.pdb") % ori_name
        if os.path.isfile(pdb_fname):
            n_found += 1
            founds[i] = pdb_fname
            shutil.copyfile(pdb_fname, os.path.join(out_dir, os.path.basename(pdb_fname)))

with open( os.path.join(out_dir, "num_to_filename.tab"), "w") as f:
    for i in sorted(founds):
        f.write(("%s\t%s\n"%(i, founds[i])))


print( n_found, i)
