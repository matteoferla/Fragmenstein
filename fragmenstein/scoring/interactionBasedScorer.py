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



class InteractionBaseScorer(_ScorerBase):


    def __init__(self, wdir, ):
        '''
        This params are generally provided through the class method computeScoreForMolecules
        :param wdir:
        '''
        super().__init__(wdir)




    def computeScoreOneMolecule(self, mol_id, mol, frags_dict, *args, **kwargs):
        '''
        :param mol_id: an id of the molecule. Only required for interoperability reasons while saving checkpoints
        :param mol: a molecule to evaluate
        :param frags_dict: a dict of fragId -> Chem.Mol to compare with bit
        :return:
        '''
        #TODO: AQUI
        return None

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
    InteractionBaseScorer.evalPipeline(initiaze_parallel_execution=True)

'''
python -m fragmenstein.scoring.xcos -i /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/compound_vs_fragments.csv -d /home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned -f /home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols  -o xcos_out.csv -p "Mpro-(\w+)\.mol" -w ~/tmp/test_dask
'''


'''
from plip.structure.preparation import PDBComplex
my_mol = PDBComplex()
my_mol.load_pdb('/home/ruben/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x3124_0A/Mpro-x3124_0A_bound.pdb')
my_mol.analyze()

my_interactions = my_mol.interaction_sets['LIG:A:1101']

my_interactions.__dict__.keys()

print([ ((interaction.resnr, interaction.reschain, interaction.restype) for interaction in my_interactions.all_hbonds_pdon]) # Prints [84, 129]

interaction.resnr_lig

'''