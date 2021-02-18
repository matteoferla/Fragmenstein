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

    mol_binary=b'\xef\xbe\xad\xde\x00\x00\x00\x00\x0c\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\r\x00\x00\x00\r\x00\x00\x00\x80\x01\x07@(\x00\x00\x00\x03\x03\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x06@(\x00\x00\x00\x03\x04\x06@h\x00\x00\x00\x03\x03\x01\x08\x00(\x00\x00\x00\x03\x02\x06\x00(\x00\x00\x00\x03\x04\x07\x00(\x00\x00\x00\x03\x03\x07\x00h\x00\x00\x00\x03\x01\x02\x07\x00(\x00\x00\x00\x03\x03\x06\x00h\x00\x00\x00\x03\x03\x01\x08\x00h\x00\x00\x00\x03\x01\x01\x0b\x00\x01h\x0c\x01\x02h\x0c\x00\x03h\x0c\x03\x04h\x0c\x02\x05h\x0c\x04\x05h\x0c\x04\x06 \x06\x07 \x07\x08.\x02\x05\x01\x00\x07\t \x08\n \n\x0b.\x02\x05\x01\x00\x0b\x0c \x14\x01\x06\x00\x03\x04\x05\x02\x01\x17\x01\x00\x00\x00\x01\x00\x00\x00\x00\rD\x8b\xda\xc17\x89!Bb\x10g\xc2\xf0\xa7\xdd\xc1\x08\xac&B\xc9\xf6g\xc2\xc3\xf5\xda\xc1\x81\x95)B\xee|l\xc2B`\xd4\xc1\xdd$\x1fB\x7f\xeaj\xc2\x93\x18\xd1\xc1V\x8e!B\xfc\xa9o\xc2\xc5 \xd4\xc1\xb4\xc8&B\xfc\xa9p\xc2T\xe3\xca\xc1u\x13\x1fB\xdd\xa4s\xc2\xbe\x9f\xc3\xc1\x17\xd9\x1aB\xd5\xf8r\xc2F\xb6\xbf\xc1j\xbc\x18B\xdd\xa4n\xc2\xa8\xc6\xc0\xc1\\\x0f\x19B5\xdew\xc2`\xe5\xc2\xc1!\xb0\x1aB9\xb4i\xc2?5\xc3\xc1}\xbf\x16Byif\xc2u\x93\xc0\xc1X\xb9\x11B\xa4\xf0g\xc2\x16'
    mol_id = 'x0020_0B_0-x0020_0B_3-x0020_0B_4-x0029_0B_2'
    mol_fname = os.path.join("/home/ruben/oxford/tools/Fragmenstein/output", mol_id, mol_id+".minimised.mol")
    from rdkit import Chem
    proposed_mols = {mol_id: (Chem.Mol(mol_binary), ["x0020_0B", "x0029_0B", "x0257_0B"] ) }
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