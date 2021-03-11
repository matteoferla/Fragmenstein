import os
from rdkit import Chem
from rdkit.Chem import PandasTools
import matplotlib.pyplot as plt
import numpy as np

sdf_fname= os.path.expanduser('~/nsp3_Bs1_x0034_0B,x0176_0B,x0183_0B,x0208_0A,x0212_0B,x0246_0B,x0276_0B,x0283_0B,x0311_0B,x0438_0B.sdf')

with open(sdf_fname, "rb") as f:

    first_line = f.readline().decode("UTF-8")
    if first_line.startswith("ver_1.2"):
        for line in f:
            if line.decode("UTF-8").startswith("$$$$"):
                break
    else:
        f.seek(0)
    df = PandasTools.LoadSDF(f, idName="molId", smilesName='smi', molColName='mol', embedProps=True, includeFingerprints=False)

print( df.shape )

props = df.columns
print(df.columns)
print(df.head())

print(df["SA_score"])

#plt.hist( np.round(df["SA_score"].values,3) ); plt.show()

print(df.shape)
df = df[ (1< df["SA_score"]< 6) & (2< df["SC_score"]< 10) ]
print(df.shape)
