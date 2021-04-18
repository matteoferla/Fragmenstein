import os
from itertools import combinations

#nsp13
# frags = "x0034_0B,x0176_0B,x0183_0B,x0208_0A,x0212_0B,x0246_0B,x0276_0B,x0283_0B,x0311_0B,x0438_0B"
# frags = "x0041_0A,x0058_0B,x0116_0B,x0276_0A,x0306_0B,x0309_0B,x0494_0B,x0499_0B"
# dataDir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"
# projectName = "nsp13-"

#dpp11
# frags = "x0032_0A,x0120_0A"
frags = "x0051_0A,x0056_0A,x0083_0A,x0087_0A,x0115_0A,x0199_0A,x0208_0A,x0228_0A,x0230_0B,x0267_0B,x0346_0A"

dataDir = "/home/ruben/oxford/myProjects/diamondOthers/fragmenstein/PgDPP11/fragalysis/aligned"
projectName="PGN_RS02895PGA-"

frags = frags.split(",")

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