from scipy.stats import gennorm

class FragmensteinScorer(): #TO be inherited by Victor?
    def __init__(self):
        pass

    def compute_ranking(self, use_old=False:

score["score_LE"] = float(score['∆∆G'] ) /(score["N_unconstrained_atoms"] + score["N_constrained_atoms"])
score["score_fragmensteinOld"] = float(score['∆∆G'] ) /5 + float(score["comRMSD"]) + score
                                                                                         ["N_unconstrained_atoms"] /5 - \
                                 score["N_constrained_atoms"] / 10


rank_weights = {'LE': 1., 'comRMSD': 2., 'atom_bonus': 2., 'novelty_penalty': 5.}
zz = (score['N_constrained_atoms'] ** 2 - 25 ** 2) / 500
atom_bonus = gennorm.pdf(zz, 5) / 0.5445622105291682
novelty_penalty = score['N_unconstrained_atoms'] / score['N_constrained_atoms']
fragmenstein_new = rank_weights['LE'] * float(score['score_LE']) + rank_weights['comRMSD'] * float(score['comRMSD']) + \
                   - rank_weights['atom_bonus'] * atom_bonus + rank_weights['novelty_penalty'] * novelty_penalty
score["score_fragmensteinNew"] = fragmenstein_new