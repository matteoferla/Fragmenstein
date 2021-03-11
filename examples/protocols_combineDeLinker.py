import numpy as np
from rdkit import Chem
from typing import Tuple, List, Dict

import itertools

from rdkit.Chem.Descriptors import HeavyAtomMolWt
from scipy.spatial import distance_matrix

from examples.protocols_combineBase import Protocol_combineBase
from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.combineMerge_DeLinkerDefault import CombineMerge_DeLinkerDefault
from fragmenstein.protocols.steps.hitsPreprocess_fragmentationBrics import HitsPreprocess_fragmentationBRICS


class Protocol_combineFragmenstein(Protocol_combineBase):

    def _check_compability(self, frag1, frag2, min_dist_thr=3,  max_dist_thr=15, min_weight=60): # Many bugs if this is 50 or less

        weight1 = HeavyAtomMolWt(frag1)
        weight2 = HeavyAtomMolWt(frag2)
        if weight1< min_weight or weight2 < min_weight:
            return False

        coords1 =  frag1.GetConformer().GetPositions()
        coords2 =  frag2.GetConformer().GetPositions()
        dist_mat = distance_matrix( coords1, coords2)

        min_val = np.min(dist_mat)
        max_val = np.max(dist_mat)

        return min_val >= min_dist_thr and max_val <= max_dist_thr

    def preprocess_fragments(self) -> Tuple[Tuple[List[Compound], Dict[str, str]]]:


        if self.preprocess_mode:
            print("preprocessing fragments")
            if self.preprocess_mode=="BRICS_decomposition":
                fragmentator = HitsPreprocess_fragmentationBRICS(self.fragments, break_mode="binary",
                                                                 random_seed=self.random_seed)
            elif self.preprocess_mode=="SMARTS_decomposition":
                raise NotImplementedError()
            else:
                raise ValueError("Error, not implemented option preprocess_mode=%s"%self.preprocess_mode)
        else:
            return None

        dict_of_frags = fragmentator.broken_fragments

        frags_keys_pairs = list( itertools.combinations( dict_of_frags.keys(), 2 ) )

        bitId_to_molId  =  fragmentator.bitId_to_molId

        #TODO: add redundancy filter
        def fragsCombin_iter():
            for key1, key2 in frags_keys_pairs:
                frags_ops1 = dict_of_frags[key1]
                frags_ops2 = dict_of_frags[key2]
                for bits1, bits2 in itertools.product(frags_ops1, frags_ops2):
                    for bit1 in bits1:
                        for bit2 in bits2:
                            if self._check_compability(bit1, bit2, min_dist_thr=5):
                                yield ( bit1, bit2 )

        # print( [ [Chem.MolToSmiles(mol) for mol in comb]  for comb in fragsCombin_iter() ] )
        # print( len(list(fragsCombin_iter())))

        return fragsCombin_iter(), bitId_to_molId

    @property
    def combiner_class(self):
        return CombineMerge_DeLinkerDefault


if __name__ == "__main__":

    parser = Protocol_combineFragmenstein.generateCmdParser(prog="protocols_combineDeLinker",
                                                            description="Combines molecules using fragmentation and DeLinker")


    args = vars(parser.parse_args())
    print(args)
    Protocol_combineFragmenstein.main(**args)
    print("\nmain DONE!\n")

    #TODO: modify fragalysis writer to accept method=DeLInker

    '''

N_CPUS=1 python -m examples.protocols_combineDeLinker -i ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -f x0176_0B x0246_0B x0438_0B -o ~/oxford/tools/Fragmenstein/output -m 10 -p BRICS_decomposition

python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json DASK_WORKER_MEMORY=4GB --ncpus 8 "/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python -m examples.protocols_combineDeLinker --n_cpus 8 -i /data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned -o /data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_nsp13_Site7_mixed -f x0116_0B x0309_0B x4094_0B"

    '''