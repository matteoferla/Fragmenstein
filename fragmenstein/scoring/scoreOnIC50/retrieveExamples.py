import os
import numpy as np
import pandas as pd

from fragmenstein.scoring.scoring_config import MPRO_RAW_DATA_DIR, RAW_SUBMISSIONS_TABLE_PATH, RAW_SUBMISSIONS_TABLE_URL

INTERESTING_FEATURES=["MW", "cLogP", "HBD", "HBA", "TPSA", "relative_solubility_at_20_uM", "relative_solubility_at_100_uM", "trypsin_IC50", "NMR_std_ratio"]
INHIBITION_FEATURES=["r_inhibition_at_20_uM", "r_inhibition_at_50_uM", "r_avg_IC50",
                     "f_inhibition_at_20_uM", "f_inhibition_at_50_uM", "f_avg_IC50"]
LABEL_FEATURE= "BINARY_LABEL"

DELTA_CONFIDENCE=0.05

def _getPosteraData():
    try:
        dataPostera = pd.read_csv(RAW_SUBMISSIONS_TABLE_PATH, error_bad_lines=False)
    except FileNotFoundError:
        dataPostera = pd.read_csv(RAW_SUBMISSIONS_TABLE_URL, error_bad_lines=False)
    return dataPostera


def retrieve_Mpro_structure_compounds(filter_out_without_structures= False, label_feature=LABEL_FEATURE):


    data = _getPosteraData()

    if filter_out_without_structures:
        data = data[(data["Structure ID"].notna())]
        compounds_with_structure = [elem.split("_")[0].split("-")[1] for elem in
                                    os.listdir(os.path.join(MPRO_RAW_DATA_DIR, "aligned"))]
        data = data[ data["Structure ID"].isin( compounds_with_structure ) ]


    #remove not assayed compounds
    data = data[ data["ASSAYED"]  ]

    data= data[ ["CID", "SMILES", 'Structure ID']+INTERESTING_FEATURES+INHIBITION_FEATURES+ ["fragments"] ]


    #remove compounds with no fragment information or known to be wrong information (x0072) and more than one inspirational hit
    data = data[ (data["fragments"].notna()) & (data["fragments"]!="x0072") & (data["fragments"].str.contains(',',case=False) ) ]


    for cname in INHIBITION_FEATURES+INTERESTING_FEATURES:  print(cname, data[cname].notna().sum())


    postive_data_mask = ((data["r_avg_IC50"]<70-DELTA_CONFIDENCE) | (data["f_avg_IC50"]<99-DELTA_CONFIDENCE) )

    data[label_feature] =  postive_data_mask
    dataF = data[ (data[label_feature].notna())]

    return dataF


if __name__=="__main__":
    data= retrieve_Mpro_structure_compounds()
    print(data) #[["Structure ID", "fragments", LABEL_FEATURE]].head(2))
    print( np.sum(data[LABEL_FEATURE]))
    print(data.shape)
'''
python -m fragmenstein.scoring.retrieveExamples
'''