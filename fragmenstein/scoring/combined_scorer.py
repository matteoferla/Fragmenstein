from typing import Dict, Tuple, List

from fragmenstein.scoring._scorer_base import _ScorerBase
from fragmenstein.scoring.propertiesScorer import PropertiesScorer


class CombineScorer(_ScorerBase):

    def __init__(self, scorers_objects_list, *args, **kwargs):
        self.scorers_objects_list = scorers_objects_list
        super().__init__(*args, **kwargs)
    @property
    def fragments_id(self):
        '''
        This is instantiation of abstract attribute
        :return:
        '''
        f_ids = []
        for scorer in self.scorers_objects_list:
            f_ids+= scorer.fragments_id

        return list( set(f_ids ) )

    def computeScoreOneMolecule(self, *args, **kwargs) -> Dict[str, Tuple[float, List[str]]]:
        '''

        :param mol_id:
        :param mol: the molecule to be scored
        :param frag_ids: the xchem ids for the inspirationl hits
        :param args:
        :param kwargs:
        :return: A dict w
        '''
        results = {}
        for scorer in self.scorers_objects_list:
            result = scorer.computeScoreOneMolecule(*args, **kwargs)
            if _ScorerBase.MOL_NAME_ID in results:
              assert results[_ScorerBase.MOL_NAME_ID] == result[_ScorerBase.MOL_NAME_ID], "Error, mismatch betweeing differnt scorers"
            results.update( result)

        return results


def test():
    import os, tempfile
    from fragmenstein.scoring.sucos import SuCOSComputer

    mol_id = 'x0020_0B_0-x0020_0B_1-x0020_0B_2-x0020_0B_3-x0020_0B_4-x0257_0B_1'
    mol_fname = os.path.join("/home/ruben/oxford/tools/Fragmenstein/output", mol_id, mol_id+".minimised.mol")
    proposed_mols = {mol_id: (mol_fname, ["x0020_0B", "x0029_0B", "x0257_0B"] ) }
    hits_root_dir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"
    fragment_id_pattern = r".*-(\w+)\.mol$"

    with tempfile.TemporaryDirectory() as tmp_score:
        su_tmp = os.path.join(tmp_score, "sucos")
        os.mkdir( su_tmp)
        scorer1 = SuCOSComputer(fragments_dir=hits_root_dir, fragment_id_pattern=fragment_id_pattern,
                                working_dir=su_tmp )
        prop_tmp = os.path.join(tmp_score, "prop")
        os.mkdir( prop_tmp)
        scorer2 = PropertiesScorer(working_dir=prop_tmp)
        scorers_list = [scorer1, scorer2]
        scores = CombineScorer.computeScoreForMolecules(proposed_mols , scorers_objects_list=scorers_list, working_dir=tmp_score)
    print(scores)


if __name__ == "__main__":
    test()

'''
python -m fragmenstein.scoring.combined_scorer
'''