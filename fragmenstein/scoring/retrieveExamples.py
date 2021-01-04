import os

import pandas as pd

from fragmenstein.scoring.scoring_config import RAW_SUBMISSIONS_TABLE_URL, MPRO_RAW_DATA_DIR

final_analysis_compounds = [ elem.split("_")[0] for elem in os.listdir( os.path.join(MPRO_RAW_DATA_DIR, "aligned") ) ]

dataPostera = pd.read_csv(RAW_SUBMISSIONS_TABLE_URL, error_bad_lines=False)
dataPostera= dataPostera[["CID", "fragments"]]
dataPostera.rename({"CID":"Compound ID"}, axis=1, inplace=True)

data = pd.read_csv(os.path.join(MPRO_RAW_DATA_DIR, 'Mpro_compound_tracker_csv.csv'), error_bad_lines=False)


print( data.head(1) )
print( data.shape )

data = data[ (data["Dataset"].notna()) ]

data = data[ data["Dataset"].isin( final_analysis_compounds ) ]

data = pd.merge( data, dataPostera, on=["Compound ID"] )

data= data[ ["SMILES", 'Dataset','Rapid Fire inhibition at 50 uM', 'Rapid Fire avg IC50 (uM)',
       'Fluorescence inhibition at 50 uM', 'Fluorescence avg IC50 (uM)',
       'NMR std ratio', 'Trypsin IC50 (uM)', "fragments"] ]


data = data[ (data["fragments"].notna()) ]


for cname in list(data.columns):
  print(cname, data[cname].isna().sum())


dataF = data[ (data["Fluorescence avg IC50 (uM)"].notna()) ]
# dataF.sort_values("Fluorescence avg IC50 (uM)")


"fragments"










'''
python -m fragmenstein.scoring.simpleScoring
'''