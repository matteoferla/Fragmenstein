import os
import sys
from typing import List, Dict, Tuple

import joblib
import pandas as pd
import pyrosetta
from joblib import Parallel
from joblib import delayed
from rdkit import Chem
from rdkit.Chem import Crippen

from fragmenstein.external import ExternalToolImporter
from fragmenstein.mpro import MProVictor
from fragmenstein.scoring.inhibitionPred.retrieveExamples import retrieve_Mpro_structure_compounds, LABEL_FEATURE
from fragmenstein.scoring.scoring_config import MPRO_RAW_DATA_DIR, FEATURES_DIR, FRAGMENSTEIN_OUT_DIR

FRAGMENSTEIN_FEATURES= ['∆∆G', '∆G_bound', '∆G_unbound', 'comRMSD']
MOLECULE_FEATURES= ['Synthetic_accessibility', 'WC_LogP' ]

class ComputeFragmenstein():

    [sascorer] = ExternalToolImporter.import_tool("DeLinker", ["sascorer"])
    def __init__(self, raw_data_path= MPRO_RAW_DATA_DIR, results_dir=FEATURES_DIR, n_jobs=1):
        '''
        :param raw_data_path: The directory where XChem data is.
        :param working_dir:
        '''

        self.raw_data_path = raw_data_path
        self.all_compounds_aligned_path = os.listdir( os.path.join(raw_data_path, "aligned") )
        self.results_dir= results_dir
        self.wdir= os.path.join(FEATURES_DIR, "working_dir")
        pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

        self.n_jobs = n_jobs


    def compute_all(self, compoundIds_to_fragmentsAndSmile: Dict[str, Tuple[List[str], str]]):
        '''
        :param compoundIds_to_fragmentsAndSmile: Dict of compound_ids-> [fragment_compound_id, smiles] to study
        :return:
        '''
        compound_id_fragments= ((key,) + compoundIds_to_fragmentsAndSmile[key] for key in compoundIds_to_fragmentsAndSmile)
        results=Parallel(n_jobs=self.n_jobs, backend="multiprocessing", batch_size=1)( delayed(self.safe_compute_stats_one_compound)(* args) for args in  compound_id_fragments )
        return results

    def safe_compute_stats_one_compound(self, compound_id:str, fragments: List[str], smiles:str):
        try:
            return self.compute_one_compound(compound_id, fragments, smiles)
        except Exception as e:
            return None

    def compute_one_compound(self, compound_id:str, fragments: List[str], smiles=None):
        '''
        :param compound_id:
        :param fragments: List of fragment codes
        :return:
        '''
        #TODO: Only works for MPro, make it more generic

        result_fname= os.path.join(self.wdir, compound_id+".joblib.pckl")
        try:
            ddict = joblib.load(result_fname)
            ddict = { ( key.decode() if hasattr(key, "decode") else key): ddict[key] for key in ddict.keys() }
            return ddict
        except (IOError, OSError):
            if not smiles:
                fullPrefix = os.path.join( self.raw_data_path, "aligned", compound_id )
                molfile = os.path.join( fullPrefix, compound_id + ".mol" )

                mol = Chem.MolFromMolFile(molfile)
                smiles = Chem.MolToSmiles(mol)
            else:
                mol = Chem.MolFromSmiles( smiles )

            MProVictor.quick_renanimation = True
            MProVictor.work_path= FRAGMENSTEIN_OUT_DIR
            victor = MProVictor.from_hit_codes(smiles= smiles,
                                               hit_codes=fragments,
                                               long_name=compound_id)
            scores = victor.summarise()
            # print( scores  )
            scores["Synthetic_accessibility"] = ComputeFragmenstein.sascorer.calculateScore(mol)
            scores["WC_LogP"] = Crippen.MolLogP( mol )
            # print(compound_id, scores)
            joblib.dump(scores, result_fname)
            victor.make_pse()
            return scores


def generateData():
    n_jobs= int( sys.argv[1])

    df = retrieve_Mpro_structure_compounds(filter_out_without_structures=False)
    id_to_y={}
    data_to_frag={}
    for i, row in df.iterrows():
        c_id = row["CID"]
        data_to_frag[ c_id ] = (row["fragments"].split(","), row["SMILES"])
        id_to_y[c_id] = row[LABEL_FEATURE]

    print("Number of compounds to compute: %s"%len(data_to_frag) )
    computer = ComputeFragmenstein(n_jobs=n_jobs)
    stats = computer.compute_all(data_to_frag)
    rows=[]

    selected_fields = FRAGMENSTEIN_FEATURES + MOLECULE_FEATURES
    for ddict in stats:
        if not ddict:
            continue
        print(ddict)
        if ddict["regarded"] ==  ['x0072']:
            continue
        try:
            rows.append(  [ddict["name"], id_to_y[ddict["name"]]] + [ ddict[key] for key in selected_fields ] )
        except KeyError:
            continue

    df= pd.DataFrame(rows, columns=["name", "label"]+selected_fields)
    print(df)
    fname_prepared_dataset = os.path.join( computer.results_dir, "features.csv")

    df.to_csv(fname_prepared_dataset, index=False)

if __name__=="__main__":
    generateData()


    '''
python -m fragmenstein.scoring.launchFragmenstein 4
    '''
