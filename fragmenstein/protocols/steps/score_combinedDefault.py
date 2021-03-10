import os
import tempfile

from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.scoring.combined_scorer import CombineScorer
from fragmenstein.scoring.interactionBasedScorer import InteractionBasedScorer
from fragmenstein.scoring.propertiesScorer import PropertiesScorer
from fragmenstein.scoring.sucos import SuCOSComputer
from fragmenstein.scoring.xcos import XcosComputer


class Score_CombinedDefault(Xchem_info):


    def __init__(self, fragments_dir, to_score_dir, fragment_id_pattern, boundPdb_id_pattern, predicted_boundPdb_id_pattern,
                 selected_fragment_ids=None, working_dir=None, *args, **kwargs):
        self.selected_fragment_ids = selected_fragment_ids
        self.to_score_dir = to_score_dir
        self.fragment_id_pattern = fragment_id_pattern
        self.fragments_dir = fragments_dir
        self.boundPdb_id_pattern = boundPdb_id_pattern
        self.predicted_boundPdb_id_pattern = predicted_boundPdb_id_pattern
        self.working_dir = working_dir

    def _getWDirScopeManager(self):
        if self.working_dir:
            class WorkingDir():
                def __enter__(self2):
                    if not os.path.isdir(self.working_dir):
                        os.mkdir(self.working_dir)
                    return self.working_dir
                def __exit__(self, exc_type, exc_val, exc_tb):
                    pass
            return WorkingDir
        else:
            return tempfile.TemporaryDirectory

    def compute_scores(self, proposed_mols, already_computed_scores=None ):
        '''

        :param proposed_mols:     proposed_mols[merge_id] = [ Chem.Mol, [fragId_1, fragId_2...] ]
        :param already_computed_scores:
        :return:
        '''

        with self._getWDirScopeManager()() as tmp:

            # scorer1 = SuCOSComputer(fragments_dir=self.fragments_dir, fragment_id_pattern=self.fragment_id_pattern, use_weights=True,
            #                         selected_fragment_ids= self.selected_fragment_ids, working_dir=tmp )
            # scorer2 = PropertiesScorer(working_dir=tmp)
            # scorer3 = XcosComputer(fragments_dir=self.fragments_dir, fragment_id_pattern=self.fragment_id_pattern,
            #                         selected_fragment_ids=self.selected_fragment_ids, working_dir=tmp )
            # scorer4 = InteractionBasedScorer(fragments_dir=self.fragments_dir, fragment_id_pattern=self.boundPdb_id_pattern,
            #                                  selected_fragment_ids=self.selected_fragment_ids, boundPdbs_to_score_dir= self.to_score_dir,
            #                                  boundPdbs_to_score_pattern= self.predicted_boundPdb_id_pattern, working_dir=tmp)

            scorers_classes = [SuCOSComputer, PropertiesScorer, XcosComputer ] #, InteractionBasedScorer]
            scorers_args = [
                dict(fragments_dir=self.fragments_dir, fragment_id_pattern=self.fragment_id_pattern, use_weights=True,
                                    selected_fragment_ids= self.selected_fragment_ids, working_dir=tmp ),
                dict(fragments_dir=self.fragments_dir, fragment_id_pattern=self.fragment_id_pattern,
                                    selected_fragment_ids=self.selected_fragment_ids, working_dir=tmp ),
                dict(fragments_dir=self.fragments_dir, fragment_id_pattern=self.fragment_id_pattern,
                                    selected_fragment_ids=self.selected_fragment_ids, working_dir=tmp ),
                dict(fragments_dir=self.fragments_dir,
                                                 fragment_id_pattern=self.boundPdb_id_pattern,
                                                 selected_fragment_ids=self.selected_fragment_ids,
                                                 boundPdbs_to_score_dir=self.to_score_dir,
                                                 boundPdbs_to_score_pattern=self.predicted_boundPdb_id_pattern,
                                                 working_dir=tmp)
            ]

            scorers_iter = list( scorer(**s_kwargs) for scorer, s_kwargs in zip(scorers_classes,scorers_args) )

            proposed_mols_dict = { mol.molId: (mol, mol.getFragIds()) for mol in proposed_mols}
            print("All scorers initialized",flush=True)
            scores_list = CombineScorer.computeScoreForMolecules(proposed_mols_dict , scorers_objects_list=scorers_iter, working_dir=tmp)


            if already_computed_scores:
                for i in range(len(scores_list)):
                    record = scores_list[i]
                    already_record = already_computed_scores[record[CombineScorer.MOL_NAME_ID]]
                    CombineScorer.update_dict(already_record, record, "_secondary" )
                    scores_list[i] = record

        # scores_dict = { record[CombineScorer.MOL_NAME_ID]: (proposed_mols[record[CombineScorer.MOL_NAME_ID]][0], record )  for record in scores_list }

            for mol, record in zip(proposed_mols, scores_list):
                mol.add_scores(record)

        return proposed_mols