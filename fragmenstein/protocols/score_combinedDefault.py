import tempfile

from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.scoring.combined_scorer import CombineScorer
from fragmenstein.scoring.interactionBasedScorer import InteractionBasedScorer
from fragmenstein.scoring.propertiesScorer import PropertiesScorer
from fragmenstein.scoring.sucos import SuCOSComputer
from fragmenstein.scoring.xcos import XcosComputer


class Score_CombinedDefault(Xchem_info):


    def __init__(self, fragments_dir, to_score_dir, fragment_id_pattern, boundPdb_id_pattern):
        self.to_score_dir = to_score_dir
        self.fragment_id_pattern = fragment_id_pattern
        self.fragments_dir = fragments_dir
        self.boundPdb_id_pattern = boundPdb_id_pattern

    def compute_scores(self, proposed_mols, already_computed_scores=None ): #TODO: add more metadata to scores
        '''

        :param proposed_mols:     proposed_mols[merge_id] = [ Chem.Mol, [fragId_1, fragId_2...] ]
        :param already_computed_scores:
        :return:
        '''
        with tempfile.TemporaryDirectory() as tmp:
            scorer1 = SuCOSComputer(fragments_dir=self.fragments_dir, fragment_id_pattern=self.fragment_id_pattern, use_weights=True, working_dir=tmp )
            scorer2 = PropertiesScorer(working_dir=tmp)
            scorer3 = XcosComputer(fragments_dir=self.fragments_dir, fragment_id_pattern=self.fragment_id_pattern, working_dir=tmp )
            scorer4 = InteractionBasedScorer( fragments_dir=self.fragments_dir, fragment_id_pattern=self.boundPdb_id_pattern,
                                              boundPdbs_to_score_dir= self.to_score_dir, boundPdbs_to_score_pattern= self.boundPdb_id_pattern, working_dir=tmp )
            scorers_list = [scorer1, scorer2, scorer3, scorer4]
            scores_list = CombineScorer.computeScoreForMolecules(proposed_mols , scorers_objects_list=scorers_list, working_dir=tmp)
            if already_computed_scores:
                for i in range(len(scores_list)):
                    record = scores_list[i]
                    already_record = already_computed_scores[record[CombineScorer.MOL_NAME_ID]]
                    CombineScorer.update_dict(already_record, record, "_secondary" )
                    scores_list[i] = already_record

        scores_dict = { record[CombineScorer.MOL_NAME_ID]: (proposed_mols[record[CombineScorer.MOL_NAME_ID]][0], record )  for record in scores_list }

        return scores_dict