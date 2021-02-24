import numpy as np
from scipy.stats import gennorm

from fragmenstein.scoring._scorer_base import _ScorerBase


class _FragmensteinScorer(): #TODO: Should this be inherited by Victor?

    old_scoring_name = _ScorerBase.SCORE_NAME_TEMPLATE%"fragmensteinOld"
    new_scoring_name = _ScorerBase.SCORE_NAME_TEMPLATE%"fragmensteinNew"
    new_rank_weights = {'LE': 1., 'comRMSD': 2., 'atom_bonus': 2., 'novelty_penalty': 5.}

    def __init__(self):
        pass

    def compute_ranking(self, use_old=False):
        record = self.summarise()
        if use_old: #TODO: should we add the computed stuff to victor __dict__?
            return _FragmensteinScorer.old_scoring_fun(record)
        else:
            return _FragmensteinScorer.new_scoring_fun(record)

    @classmethod
    def old_scoring_fun(cls, record):
        new_record = record.copy()
        new_record[cls.old_scoring_name] = float(record['∆∆G'] ) /5 + float(record["comRMSD"]) + \
                                           record["N_unconstrained_atoms"] /5 - record["N_constrained_atoms"] / 10
        return _FragmensteinScorer.old_scoring_name, new_record


    @classmethod
    def new_scoring_fun(cls, record):
        new_record = record.copy()
        new_record[_ScorerBase.SCORE_NAME_TEMPLATE%'LE'] = float(record['∆∆G'] ) /(record["N_unconstrained_atoms"] + record["N_constrained_atoms"])
        new_record[_ScorerBase.SCORE_NAME_TEMPLATE%'comRMSD'] = new_record['comRMSD']
        zz = (new_record['N_constrained_atoms'] ** 2 - 25 ** 2) / 500
        atom_bonus = gennorm.pdf(zz, 5) / 0.5445622105291682

        novelty_penalty = (1+new_record['N_unconstrained_atoms'])/ (1+new_record['N_constrained_atoms'])

        fragmenstein_new = cls.new_rank_weights['LE'] * float(new_record['LE_score']) + cls.new_rank_weights['comRMSD'] * float(new_record['comRMSD']) + \
                           - cls.new_rank_weights['atom_bonus'] * atom_bonus + cls.new_rank_weights['novelty_penalty'] * novelty_penalty

        new_record[_FragmensteinScorer.new_scoring_name] = fragmenstein_new
        return _FragmensteinScorer.new_scoring_name, new_record
