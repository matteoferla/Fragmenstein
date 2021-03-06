#!/usr/bin/env python

"""
This is a reimplementation/modification of the XCos code that was originally written by Warren Thompson
<warren.thompson@diamond.ac.uk> taken from (https://github.com/InformaticsMatters/pipelines/blob/master/src/python/pipelines/xchem/xcos.py)
"""
import os
import re

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS

from fragmenstein.protocols.steps.hitsPreprocess_fragmentationBrics import split_mol_to_brics_bits
from fragmenstein.scoring._scorer_base import _ScorerBase
from fragmenstein.scoring.cos_like_base import _COSLikeBase


class XcosComputer(_COSLikeBase): #TODO: Refactor xcos tu use cos_like_base


    def __init__(self, do_fragments_of_fragments=False, perBit_score_threshold=None, perFragment_score_threshold=None, *args, **kwargs):
        '''
        This params are generally provided through the class method computeScoreForMolecules directly obtained from cmd parser
        '''


        super().__init__( *args, **kwargs)

        self.do_fragments_of_fragments = do_fragments_of_fragments
        if do_fragments_of_fragments:

            self.fragments_dict = self.fragment_fragments( self.fragments_dict )

            perBit_score_threshold = perBit_score_threshold if perBit_score_threshold else 0.2
            perFragment_score_threshold = perFragment_score_threshold if perFragment_score_threshold else 0.5
        else:
            perBit_score_threshold = perBit_score_threshold if perBit_score_threshold else 0.4
            perFragment_score_threshold = perFragment_score_threshold if perFragment_score_threshold else 2.0

        self.perBit_score_threshold = perBit_score_threshold
        self.perFragment_score_threshold = perFragment_score_threshold




    def fragment_fragments(self, fragments=None):

        if fragments is None:
            fragments = self.fragments_dict

        new_fragments_dict= {}
        for frag_id, fragment in fragments.items():
            for i, bit in  enumerate(split_mol_to_brics_bits(fragment)):
                bit = Chem.DeleteSubstructs(bit, Chem.MolFromSmarts('[#0]'))
                new_fragments_dict[frag_id+"_%d"%i] = bit

        return new_fragments_dict

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
            protrude_score = self.getShapeScore( bit, frag_mol, vdwScale=0.2)
            # print(Chem.MolToSmiles(bit), Chem.MolToSmiles(frag_mol), protrude_score)

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

            if self.do_fragments_of_fragments:
                modifyFragId = lambda x : self._getFragIdFromSubFrag(x)
            else:
                modifyFragId = lambda  x: x
            scores_iter = map(lambda fragId_frag: ( modifyFragId(fragId_frag[0]), self.computeScoreBitFragPair(bit, n_bitAtoms, n_bitAtoms_without_wildCard, fragId_frag[1])),
                              frags_dict.items())

            scores_iter = filter(lambda x: x[-1] > self.perFragment_score_threshold, scores_iter)
            scores_list = list(scores_iter  )
            if len(scores_list)==0:
                return None
            else:
                return scores_list
        else:
            return None

    def _getFragIdFromSubFrag(self, subfragment_id):
        return  re.sub("_\d+$","",subfragment_id)

    def computeScoreOneMolecule(self, mol_id, mol, frag_ids, *args, **kwargs):
        '''
        :param mol_id: an id of the molecule.
        :param mol: a molecule to evaluate
        :return:
        '''

        # Get the bits
        compound_bits = split_mol_to_brics_bits(mol)

        if self.do_fragments_of_fragments:
            def checkInDict(key, query_dict):
                return self._getFragIdFromSubFrag(key) in query_dict
        else:
            def checkInDict(key, query_dict):
                return key in query_dict

        current_frag_dict = {key: self.fragments_dict[key] for key in self.fragments_dict if checkInDict(key,frag_ids) }
        fragName_score_per_bit = map(lambda bit: self.computeScoreOneBit(bit, current_frag_dict), compound_bits)
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

        partial_results = {_ScorerBase.MOL_NAME_ID: mol_id, _ScorerBase.SCORE_NAME_TEMPLATE%"xcos": score, _ScorerBase.FRAGMENTS_ID: list(fragments) }
        return partial_results

    @classmethod
    def parseCmd(cls):
        description = "XCos scoring with RDKit"
        additional_args = [
                            ('-f', '--fragments_dir', dict(required=True, type=str, help='Directory with a file for each fragment to '
                                                                                             'compare. Mol format required. Fragments can also be located in subdirectories)')),

                            ('-t', '--perBit_score_threshold', dict(type=float, default=0.4,
                                                                   help='Minimum shape overlay and feature map score required for scoring a bit to a fragment. Default: %(default)s')),

                            ('-g', '--perFragment_score_threshold', dict(type=float, default=2.0,
                                                                        help='Minimun per bit score summation to consider a hit as succesful'
                                                                             '. Default: %(default)s'))
                           ]
        return _ScorerBase.parseCmd(description, additional_args)




def test():
    #TODO: Remove this
    mol_id = 'x0020_0B_0-x0020_0B_1-x0020_0B_2-x0020_0B_3-x0020_0B_4-x0257_0B_1'
    mol_fname = os.path.join("/home/ruben/oxford/tools/Fragmenstein/output", mol_id, mol_id+".minimised.mol")
    proposed_mols = {"mol_id": (mol_fname, ["x0020_0B", "x0029_0B", "x0257_0B"] ) }
    hits_root_dir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"
    fragment_id_pattern = r".*-(\w+)\.mol$"
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_score:
        scores = XcosComputer.computeScoreForMolecules(proposed_mols, fragments_dir=hits_root_dir, do_fragments_of_fragments=True,
                                                       fragment_id_pattern=fragment_id_pattern, working_dir=tmp_score)
    print(scores)

if __name__ == "__main__":
    # import sys; test(); sys.exit(0)
    results = XcosComputer.evalPipeline(initiaze_parallel_execution=True)
    print("RESULTS")
    print(results)

'''
python -m fragmenstein.scoring.xcos -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/hit_mols  -o compound-set_xcosRSG.csv -s  compound-set_xcosRSG.sdf -p "Mpro-(\w+)\.mol" -w ~/tmp/test_dask

'''