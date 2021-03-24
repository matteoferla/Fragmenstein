import os

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
import matplotlib.pyplot as plt
import numpy as np

from fragmenstein.external.uploadToFragalysis.fragalysisFormater import FragalysisFormater
from fragmenstein.scoring.scorer_labels import checkIfNameIsScore, SCORE_NAME_TEMPLATE

# sdf_fname= os.path.expanduser('~/nsp13_Bs1_x0034_0B,x0176_0B,x0183_0B,x0208_0A,x0212_0B,x0246_0B,x0276_0B,x0283_0B,x0311_0B,x0438_0B.sdf')
# # sdf_fname= os.path.expanduser('~/nsp3_Bs1_x0041_0A,x0058_0B,x0176_0A,x0309_0A,x0494_0B.sdf')
#
# sdf_filtered_fname = os.path.expanduser('~/nsp13_Bs1_filtered.sdf')


sdf_fname= os.path.expanduser(sys.argv[1])

sdf_filtered_fname = None

if len( sys.argv) ==3:
    sdf_filtered_fname = os.path.expanduser(sys.argv[2])

with open(sdf_fname, "rb") as f:

    first_line = f.readline().decode("UTF-8")
    if first_line.startswith("ver_1.2"):
        for line in f:
            if line.decode("UTF-8").startswith("$$$$"):
                break
    else:
        f.seek(0)
    df = PandasTools.LoadSDF(f, idName="molId", smilesName='smi', molColName='mol', embedProps=True, includeFingerprints=False)

df[SCORE_NAME_TEMPLATE%"molMass"] = df["mol"].map( Descriptors.MolWt)


print( df.shape )
props = df.columns

for prop in props:
    if checkIfNameIsScore(prop):
       df[prop] = df[prop].astype(np.float32)
       # print(prop); plt.title(prop); plt.hist( df[prop], label= prop ); plt.show()
       # print(prop); plt.boxplot( df[prop] ); plt.show()

print(df.columns)

print(df.shape)


df = df.query( " comRMSD_score < 1.")
print("rmsd", df.shape)

df = df.query( " 100 < molMass_score < 500")
print("mass", df.shape)

# prop="SA_score"; print(prop); plt.title(prop); plt.hist( df[prop], label= prop ); plt.show()

df = df.query( " 0.1 < SA_score < 7 and 0.2 < SC_score < 3.5")
print("Synthetic accesibility", df.shape)

df = df.query( " -1 < aLogP_score < 4 ")
print("log_p", df.shape)

df = df.query( " polarSurfaceArea_score <=140 ")
print("polar surface", df.shape)

df = df.query( " 2 <rotableBonds_score <= 10 ")
print("rotable bonds", df.shape)

df = df.query( " meanSuCosW_score > 0.2 ")
print("suCos", df.shape)

df = df.query( " xcos_score >0.1 ")
print("xcos", df.shape)

df = df.query( " fragmensteinOld_score < 5 ")
print("fragmesteinOld", df.shape)

df = df.query( " fragmensteinNew_score < 130 ")
print("fragmesteinNew", df.shape)

# df = df.query( " plipGlobalPreser_score > 0.1 ")
# print("plip", df.shape)

# for mol in df["mol"]:
#     print( mol.GetPropsAsDict() )
#     input(mol)

# df = df.sort_values(by="rotableBonds_score").iloc[:40]

if sdf_filtered_fname:
    FragalysisFormater().write_molsList_to_sdf(sdf_filtered_fname,  df["mol"])

'''
python -m examples.analyze in.sdf out.sdf
'''
