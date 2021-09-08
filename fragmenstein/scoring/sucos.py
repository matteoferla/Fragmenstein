#!/usr/bin/env python
import json
import os
import numpy as np

from fragmenstein.scoring._scorer_base import _ScorerBase
from fragmenstein.scoring.cos_like_base import _COSLikeBase
from fragmenstein.utils.config_manager import ConfigManager


class SuCOSComputer(_COSLikeBase):


    def __init__(self, frags_thr=0.1, use_weights=False, ncpus = int(ConfigManager.N_CPUS), *args, **kwargs):
        '''
        This params are generally provided through the class method computeScoreForMolecules directly obtained from cmd parser
        '''
        super().__init__(*args, **kwargs)
        self.frags_thr = frags_thr
        self.use_weights = use_weights
        if use_weights:
            if isinstance(use_weights, dict):
                self.weights = use_weights
            elif isinstance(use_weights, str) and os.path.exists(use_weights):
                with open(use_weights) as f:
                    self.weights = json.load(f)
            elif isinstance(use_weights, bool):
                print("Compute occupancy weights")
                frag_ids, frags = zip(* self.fragments_dict.items() )
                weights = self.compute_occupancy_weights( frags, n_cpus = ncpus )
                self.weights = { frag_id: w for frag_id, w in zip(frag_ids, weights) }
                print("Occupancy weights computed")
            else:
                raise ValueError("Error, bad option for use_weights. Found (%s)"%str(use_weights))


    def computeScoreOneMolecule(self, mol_id, mol, frag_ids, *args, **kwargs):
        '''
        :param mol_id: an id of the molecule.
        :param mol: a molecule to evaluate
        :return:
        '''

        # Get the bits


        current_fragments = [ self.fragments_dict[key] for key in frag_ids ]
        per_fragment_score = list(map(lambda x: self.computeSuCOS(mol, x), current_fragments))
        if self.use_weights:
            per_fragment_score = [ score* self.weights[frag_id] for score, frag_id in zip(per_fragment_score, frag_ids) if score > self.frags_thr]
            sucosName = "SuCosW"
        else:
            sucosName="SuCos"

        selected_frags = [ frag_id for score, frag_id in zip(per_fragment_score, frag_ids) if score > self.frags_thr]
        if len(selected_frags) >0:
            mean_score = np.mean(per_fragment_score)
            max_score = np.max(per_fragment_score)
            min_score = np.min(per_fragment_score)
            sum_sucos = np.sum(per_fragment_score)
        else: #in case no fragments, just assign nan
            mean_score = np.nan
            max_score = np.nan
            min_score = np.nan
            sum_sucos = np.nan

        partial_results = {_ScorerBase.MOL_NAME_ID: mol_id,
                           _ScorerBase.SCORE_NAME_TEMPLATE%("mean"+sucosName): mean_score,
                           _ScorerBase.SCORE_NAME_TEMPLATE%("max"+sucosName): max_score,
                           _ScorerBase.SCORE_NAME_TEMPLATE%("min"+sucosName): min_score,
                           _ScorerBase.SCORE_NAME_TEMPLATE%("sum" + sucosName): sum_sucos,
                           _ScorerBase.FRAGMENTS_ID: sorted(selected_frags)}
        return partial_results

    @classmethod
    def parseCmd(cls):
        description = "SuCOS scoring of fragments"
        additional_args = [
                            ('-f', '--fragments_dir', dict(required=True, type=str, help='Directory with a file for each fragment to '
                                                                                             'compare. Mol format required. Fragments can also be located in subdirectories)')),
                            ( '--use_weights', dict(required=False, action="store_true", help='Use a weighted version based on fragment occupancy fraction')),

                           ]
        return super().parseCmd(description, additional_args)

def test():
    from rdkit import Chem
    mol_bin = b'\xef\xbe\xad\xde\x00\x00\x00\x00\x0c\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x13\x00\x00\x00\x13\x00\x00\x00\x80\x01\x06\x00`\x00\x00\x00\x02\x02\x06\x00`\x00\x00\x00\x02\x02\x07\x00`\x00\x00\x00\x02\x01\x06@(\x00\x00\x00\x03\x04\x06@(\x00\x00\x00\x03\x04\x06@h\x00\x00\x00\x03\x03\x01\t\x00 \x00\x00\x00\x01\x10\x00 \x00\x00\x00\x06\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x06\x00`\x00\x00\x00\x01\x03\x08\x00(\x00\x00\x00\x03\x02\x08\x00(\x00\x00\x00\x03\x02\x06@(\x00\x00\x00\x03\x04\x06\x00`\x00\x00\x00\x02\x02\x10\x00 \x00\x00\x00\x06\x08\x00(\x00\x00\x00\x03\x02\x08\x00(\x00\x00\x00\x03\x02\x07\x00`\x00\x00\x00\x01\x02\x0b\x00\x01\x00\x01\x02\x00\x00\x03\x00\x03\x04h\x0c\x03\x05h\x0c\x04\x06\x00\x07\x02\x00\x04\x08h\x0c\x05\th\x0c\x07\n\x00\x07\x0b\x08\x02\x07\x0c\x08\x02\x08\rh\x0c\t\rh\x0c\r\x0e\x00\x0e\x0f\x00\x0f\x10\x08\x02\x0f\x11\x08\x02\x0f\x12\x00\x14\x01\x06\x03\x05\t\r\x08\x04\x17\x01\x00\x00\x00\x01\x00\x00\x00\x00\x13\xa0\x1aE\xc1\xac\x1c\xaa@\x04\x96\x9a\xc2\x04VT\xc1\xecQ\xb4@\xc1J\x98\xc2J\x0cL\xc1}?\xd9@V\xce\x96\xc2\xdfO3\xc1\x83\xc0\x8a@\xc3\xf5\x99\xc2;\xdf3\xc1\xe9&Y@N\xe2\x97\xc2}?!\xc1w\xbe\x8b@J\x8c\x9b\xc2\xe5\xd0D\xc1\x08\xacT@\x0cB\x96\xc27\x89G\xc1\x85\xeb\xd1@J\x8c\x93\xc2\x91\xed"\xc1b\x10 @\xcfw\x97\xc2X9\x10\xc1)\\_@`%\x9b\xc2\xf4\xfd\\\xc1\xd7\xa3\xf0@\xb02\x92\xc2\x8bl3\xc1\xfa~\xe6@\x8b\xec\x92\xc2\xc3\xf5J\xc1\xb8\x1e\xa5@1\x08\x93\xc2\xfe\xd4\x10\xc15^"@\xc5 \x99\xc2\xc3\xf5\xfc\xc0#\xdb\xc9?\xc7\xcb\x98\xc2q=\x08\xc1\x8f\xc2\xf5\xbdF\xf6\x98\xc2\x1b/\x0f\xc1\xe5\xd0\xe2\xbeN\xa2\x9b\xc2\xe7\xfb\x15\xc17\x89\xc1\xbe\xb0\xb2\x96\xc2\x1b/\xe5\xc0\x87\x16\x89\xbf?u\x98\xc2\x16'
    proposed_mols = {"mol_id": (Chem.Mol(mol_bin), ["x0183_0B", "x0276_0B"] ) }
    hits_root_dir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"
    fragment_id_pattern = r".*-(\w+)\.mol$"
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_score:
        scores = SuCOSComputer.computeScoreForMolecules(proposed_mols, fragments_dir=hits_root_dir, use_weights=False,
                                                       fragment_id_pattern=fragment_id_pattern, working_dir=tmp_score)
    print(scores)
    print("test_done")

if __name__ == "__main__":
    # import sys; test(); sys.exit(0)

    results = SuCOSComputer.evalPipeline(initiaze_parallel_execution=True)
    print(results)

'''
python -m fragmenstein.scoring.sucos -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/hit_mols_full_name  -o compound-set_xcosRSG.csv -s  compound-set_xcosRSG.sdf -p "Mpro-(\w+)\.mol" -w ~/tmp/test_dask

N_CPUS=2 python -m fragmenstein.scoring.sucos -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/hit_mols_full_name  -o compound-set_xcosRSG.csv -s  compound-set_xcosRSG.sdf -p "Mpro-(\w+)\.mol" -w ~/tmp/test_dask --use_weights


'''