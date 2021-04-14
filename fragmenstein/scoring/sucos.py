#!/usr/bin/env python
import json
import os
import numpy as np

from fragmenstein.scoring._scorer_base import _ScorerBase
from fragmenstein.scoring.cos_like_base import _COSLikeBase



class SuCOSComputer(_COSLikeBase):


    def __init__(self, frags_thr=0.1, use_weights=False,  *args, **kwargs):
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
                frag_ids, frags = zip(* self.fragments_dict.items() )
                weights = self.compute_occupancy_weights( frags, )
                self.weights = { frag_id: w for frag_id, w in zip(frag_ids, weights) }
            else:
                raise ValueError("Error, bad option for use_weights. Found (%s)"%str(use_weights))

    def computeSuCOS(self, mol, frag):
        return 0.5*self.getFeatureMapScore(frag, mol) + 0.5*self.getShapeScore(frag, mol)

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

        try:
            mean_score = np.mean(per_fragment_score)
            max_score = np.max(per_fragment_score)
            min_score = np.min(per_fragment_score)
            sum_sucos = np.sum(per_fragment_score)
        except ValueError: #in case no fragments, just assign nan
            mean_score = np.nan
            max_score = np.nan
            min_score = np.nan
            sum_sucos = np.sum

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
    mol_id = 'x0020_0B_0-x0020_0B_3-x0029_0B_0-x0029_0B_1-x0029_0B_2-x0257_0B_0-x0257_0B_1'
    mol_fname = os.path.join("/home/ruben/oxford/tools/Fragmenstein/output", mol_id, mol_id+".minimised.mol")
    proposed_mols = {"mol_id": (mol_fname, ["x0020_0B", "x0029_0B", "x0257_0B"] ) }
    hits_root_dir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"
    fragment_id_pattern = r".*-(\w+)\.mol$"
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_score:
        scores = SuCOSComputer.computeScoreForMolecules(proposed_mols, fragments_dir=hits_root_dir, use_weights=False,
                                                       fragment_id_pattern=fragment_id_pattern, working_dir=tmp_score)
    print(scores)

if __name__ == "__main__":
    # import sys; test(); sys.exit(0)

    results = SuCOSComputer.evalPipeline(initiaze_parallel_execution=True)
    print(results)

'''
python -m fragmenstein.scoring.sucos -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/hit_mols_full_name  -o compound-set_xcosRSG.csv -s  compound-set_xcosRSG.sdf -p "Mpro-(\w+)\.mol" -w ~/tmp/test_dask

N_CPUS=2 python -m fragmenstein.scoring.sucos -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/hit_mols_full_name  -o compound-set_xcosRSG.csv -s  compound-set_xcosRSG.sdf -p "Mpro-(\w+)\.mol" -w ~/tmp/test_dask --use_weights


'''