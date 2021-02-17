#!/usr/bin/env python

"""
Scoring based on properties such as synthetic accessibility
"""

from rdkit import Chem

from rdkit.Chem.Descriptors import ExactMolWt

from fragmenstein.external import ExternalToolImporter
from fragmenstein.external.scscore.SCScoreWrapper import SCScoreWrapper
from fragmenstein.scoring._scorer_base import _ScorerBase


class PropertiesScorer(_ScorerBase):


    def __init__(self, *args, **kwargs):
        '''
        This params are generally provided through the class method computeScoreForMolecules
        args/kwargs must contain working_directory

        '''

        [sascorer] = ExternalToolImporter.import_tool("DeLinker", ["sascorer"])
        self.sascorer = sascorer

        self.scsw = SCScoreWrapper()

        super().__init__( *args, **kwargs)

    @property
    def fragments_id(self):
        '''
        This is instantiation of abstract attribute
        :return: the list of all fragment ids that have been loaded
        '''
        return None

    def computeScoreOneMolecule(self, mol_id, mol, frag_ids, *args, **kwargs):
        '''
        :param mol_id: an id of the molecule.
        :param mol: a molecule to evaluate
        :param frag_ids: ignored. Included for compatibility reasons
        :return:
        '''
        sascore = self.sascorer.calculateScore(mol)
        scscore =  self.scsw.compute_score(Chem.MolToSmiles(mol))[0]
        descriptors =  Chem.QED.properties(mol)
        partial_results = {_ScorerBase.MOL_NAME_ID: mol_id, _ScorerBase.SCORE_NAME_TEMPLATE%"SA": sascore,
                           _ScorerBase.SCORE_NAME_TEMPLATE%"aLogP":descriptors.ALOGP,  _ScorerBase.SCORE_NAME_TEMPLATE%"hbondD":descriptors.HBD,
                           _ScorerBase.SCORE_NAME_TEMPLATE%"hbondA":descriptors.HBA, _ScorerBase.SCORE_NAME_TEMPLATE%"polarSurfaceArea":descriptors.PSA,
                           _ScorerBase.SCORE_NAME_TEMPLATE%"rotableBonds":descriptors.ROTB, _ScorerBase.SCORE_NAME_TEMPLATE%"SC":scscore}
        return partial_results

    @classmethod
    def parseCmd(cls):
        description = "Computes drug-related properties for molecules"
        return _ScorerBase.parseCmd(description)


def test():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_score:
        ps = PropertiesScorer(working_dir=tmp_score)
        mol = Chem.MolFromSmiles("CC")
        results = ps.computeScoreOneMolecule("one_id", mol, None)
        print( results )

if __name__ == "__main__":
    # test() ; input("done")
    results = PropertiesScorer.evalPipeline(initiaze_parallel_execution=True)
    print(results)

'''
python -m fragmenstein.scoring.propertiesScorer -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/hit_mols  -o compound-set_properties.csv -s  compound-set_properties.sdf  -w ~/tmp/test_dask
'''
