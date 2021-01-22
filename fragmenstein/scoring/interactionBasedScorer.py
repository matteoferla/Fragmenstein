#!/usr/bin/env python

"""
This is a reimplementation/modification of the XCos code that was originally written by Warren Thompson
<warren.thompson@diamond.ac.uk> taken from (https://github.com/InformaticsMatters/pipelines/blob/master/src/python/pipelines/xchem/xcos.py)
"""
import argparse
import os

from rdkit import Chem
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit import RDConfig
from rdkit.Chem import rdFMCS

import numpy as np

from fragmenstein.scoring._scorer_base import _ScorerBase



class XcosComputer(_ScorerBase):


    def __init__(self, wdir, perBit_score_threshold=0.4, perFragment_score_threshold=6.0, feature_factory = None):
        '''
        This params are generally provided through the class method computeScoreForMolecules
        :param wdir:
        :param perBit_score_threshold:
        :param perFragment_score_threshold:
        :param feature_factory: AllChem.BuildFeatureFactor. Should be created within subprocess for parallelism.
        '''
        super().__init__(wdir)
        self.perBit_score_threshold = perBit_score_threshold
        self.perFragment_score_threshold = perFragment_score_threshold
        if feature_factory is None:
            feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
        self.feature_factory = feature_factory
        self.fmParams = {k: FeatMaps.FeatMapParams() for k in feature_factory.GetFeatureFamilies()}
        self.keep_featnames = list( self.fmParams.keys() )

    def getFeatureMapScore(self, small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.All):
        try:
            featLists = []
            for m in [small_m, large_m]:
                rawFeats = self.feature_factory.GetFeaturesForMol(m)
                # filter that list down to only include the ones we're intereted in
                featLists.append([f for f in rawFeats if f.GetFamily() in self.keep_featnames])
            fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params= self.fmParams) for x in featLists]
            fms[0].scoreMode = score_mode
            fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
            return fm_score
        except ZeroDivisionError:
            return 0

    def computeScoreBitFragPair(self, bit, n_bitAtoms, n_bitAtoms_without_wildCard, frag_mol):
        '''

        :param bit: Piece of the molecule to score
        :param n_bitAtoms: number of atoms of the bit (provided for performance reasons )
        :param n_bitAtoms_without_wildCard: number of atoms of the bit excluding wildCard atoms (provided for performance reasons )
        :param frag_mol: a fragment to evaluate against the bit
        :return: Tuple[ frag_name, final_score ]
        '''

        # Get frag name for linking to score

        # Score only if some common structure shared between bit and fragment.
        # Check if MCS yield > 0 atoms
        mcs_match = rdFMCS.FindMCS([bit, frag_mol], ringMatchesRingOnly=True, matchValences=True)

        # Get mcs_mol from mcs_match
        mcs_mol = Chem.MolFromSmarts(mcs_match.smartsString)

        # check if frag has MCS mol
        mcs_test = frag_mol.HasSubstructMatch(mcs_mol)

        if mcs_test:
            # Change van der Waals radius scale for stricter overlay
            protrude_dist = rdShapeHelpers.ShapeProtrudeDist(bit, frag_mol, allowReordering=False, vdwScale=0.2)
            protrude_dist = np.clip(protrude_dist, 0, 1)

            protrude_score = 1 - protrude_dist

            if protrude_score > self.perBit_score_threshold:
                fm_score = self.getFeatureMapScore(bit, frag_mol)
                fm_score = np.clip(fm_score, 0, 1)
                if fm_score < self.perBit_score_threshold:
                    fm_score = 0
            else:
                fm_score = 0
        else:
            protrude_score = 0
            fm_score = 0

        final_score =  0.5 * (fm_score * n_bitAtoms_without_wildCard) + 0.5 * (protrude_score * n_bitAtoms)

        return final_score


    def computeScoreOneBit(self, bit, frags_dict):
        '''
        compute the score associated to input molecule bit and select mathcing fragments
        :param bit: A piece of the molecule to evaluate
        :param frags_dict: a dict of frag_id -> Chem.Mol to compare with bit
        :return: an iterator  ( (frag_id , score ) or None  if no suitable bit provided
        '''
        # Let's remove wildcard atoms
        # Removing wildcard atoms does not impact feat score but does lower shape overlay
        # For scoring should multiply feat score by number of non-wilcard atoms and use
        # all atoms including wildcard for shape overlay
        bit_without_wildCard_atoms = Chem.DeleteSubstructs(bit, Chem.MolFromSmarts('[#0]'))

        # Let's only score bits that have more than one atom (do not count wildcard atoms)
        # Get number of bit atoms without wildcard atoms
        n_bitAtoms_without_wildCard = bit_without_wildCard_atoms.GetNumAtoms()

        # Get number of bit atoms
        n_bitAtoms = bit.GetNumAtoms()

        # Only score if enough info in bit to describe a vector - this will bias against
        # cases where frag has long aliphatic chain
        if n_bitAtoms_without_wildCard > 1:
            scores_iter = map(lambda fragId_frag: (fragId_frag[0], self.computeScoreBitFragPair(bit, n_bitAtoms, n_bitAtoms_without_wildCard, fragId_frag[1])),
                              frags_dict.items())

            scores_iter = filter(lambda x: x[-1] > self.perFragment_score_threshold, scores_iter)
            scores_list = list(scores_iter  )
            if len(scores_list)==0:
                return None
            else:
                return scores_list
        else:
            return None


    def computeScoreOneMolecule(self, mol_id, mol, frags_dict, *args, **kwargs):
        '''
        :param mol_id: an id of the molecule. Only required for interoperability reasons while saving checkpoints
        :param mol: a molecule to evaluate
        :param frags_dict: a dict of fragId -> Chem.Mol to compare with bit
        :return:
        '''
        # Get the bits
        compound_bits = self.splitMolToBits(mol)
        fragName_score_per_bit = map(lambda bit: self.computeScoreOneBit(bit, frags_dict), compound_bits)

        fragments = set([])
        score = 0
        for bit_result in fragName_score_per_bit:
            if bit_result is not None:
                frag_names, scores = zip(*bit_result)
                fragments = fragments.union( frag_names )
                score += max(scores)
            else:
                #TODO: Penalty
                pass

        partial_results = {"mol_name": mol_id, "score": score, "fragments": list(fragments) }
        return partial_results

    @classmethod
    def parseCmd(cls):
        description = "XCos scoring with RDKit"
        additional_args = [('-t', '--perBit_score_threshold', dict(type=float, default=0.4,
                                                                   help='Minimum shape overlay and feature map score required for scoring a bit to a fragment. Default: %(default)s')),

                           ('-g', '--perFragment_score_threshold', dict(type=float, default=2.0,
                                                                        help='Minimun per bit score summation to consider a hit as succesful'
                                                                             '. Default: %(default)s'))
                           ]
        return _ScorerBase.parseCmd(description, additional_args)


if __name__ == "__main__":
    XcosComputer.evalPipeline(initiaze_parallel_execution=True)

'''
python -m fragmenstein.scoring.xcos -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols  -o xcos_out.csv -p "Mpro-(\w+)\.mol" -w ~/tmp/test_dask
'''