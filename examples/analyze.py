import os

import sys
import argparse
from zipfile import ZipFile

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
import matplotlib.pyplot as plt
import numpy as np

from fragmenstein.external.uploadToFragalysis.fragalysisFormater import FragalysisFormater
from fragmenstein.scoring.scorer_labels import checkIfNameIsScore, SCORE_NAME_TEMPLATE

parser = argparse.ArgumentParser("analyze_sdf")
parser.add_argument("-i", "--input", type=str, required=True, help="input sdf file")
parser.add_argument("-o", "--output", type=str, required=False, default=None, help="output sdf file")
parser.add_argument("-p", "--pdbZipIn", type=str, required=False, default=None, help="input pdb zip file")
parser.add_argument("-b", "--pdbZipOut", type=str, required=False, default=None, help="output pdb zip file")

args = parser.parse_args()
sdf_fname= os.path.expanduser(args.input)

sdf_filtered_fname = os.path.expanduser(args.output)  if args.output else None
inPdbs_fname = os.path.expanduser(args.pdbZipIn)  if args.pdbZipIn else None
outPdbs_fname = os.path.expanduser(args.pdbZipOut)  if args.pdbZipOut else None


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



# df = df.sort_values(by="rotableBonds_score").iloc[:40]

if sdf_filtered_fname:
    FragalysisFormater().write_molsList_to_sdf(sdf_filtered_fname,  df["mol"])


#TODO: if zip file provided extract the associated molecules https://docs.python.org/3/library/zipfile.html#zipfile.ZipFile.open

if inPdbs_fname:
    assert outPdbs_fname is not None, "Error if inPdbs_fname is not None, outPdbs_fname should be provided"
    selectedIds = set( [ mol.GetProp("original_name") for mol in df["mol"]] )
    found_pdbs = set([])
    summary_fname = None
    with ZipFile( inPdbs_fname ) as f_in, ZipFile( outPdbs_fname, "w" ) as f_out:
        for fname in f_in.namelist():
            if fname.endswith(".pdb"):
                basename = os.path.basename(fname)
                molId = basename.split(".")[0]
                if molId in selectedIds:
                    found_pdbs.add(basename)
                    file_str = f_in.read( fname )
                    f_out.writestr(fname, file_str)
            elif fname.endswith(".txt"):
                summary_fname = fname
        found_lines = []
        for line in  f_in.read( summary_fname).decode("utf-8").split("\n"):
            molName, pdbFname = line.split()
            if pdbFname in found_pdbs:
                found_lines.append( line )
        f_out.writestr(summary_fname, "\n".join(found_lines))

'''
python -m examples.analyze  -i in.sdf -o out.sdf
'''
