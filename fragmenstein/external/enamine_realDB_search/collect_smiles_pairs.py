import os
from itertools import combinations

frags = "x0034_0B,x0176_0B,x0183_0B,x0208_0A,x0212_0B,x0246_0B,x0276_0B,x0283_0B,x0311_0B,x0438_0B"
frags = frags.split(",")
projectName = "nsp13-"
dataDir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"

smiles_list = set([])
for fid in frags:
  smiFname = os.path.join(dataDir, projectName+fid, projectName+fid+"_smiles.txt")
  with open(smiFname) as f:
    smi = f.read().strip()
    smiles_list.add( smi )

smiles_list = sorted(smiles_list)

# print(dict(zip(frags, smiles_list)))
# print("#%d"%len(smiles_list))
for x,y in combinations(smiles_list, 2):
  print("%s,%s"%(x,y))