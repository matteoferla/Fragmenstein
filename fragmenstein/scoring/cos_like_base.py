#!/usr/bin/env python

"""
This is a base clase for the COS-like scores such as xcos, suCOS, etc. Those scores are based on shape and chemical complementarity.
"""

import os

import mrcfile

import numpy as np

from rdkit.Geometry import rdGeometry

from fragmenstein.external import ExternalToolImporter
from fragmenstein.scoring._scorer_base import _ScorerBase
from fragmenstein.utils.io_utils import load_mol, load_files_as_mols



sucos_class = ExternalToolImporter.import_tool("sucos", ["sucos"])[0].SuCOS

class _COSLikeBase(_ScorerBase, sucos_class):
    """
    This is a base clase for the COS-like scores such as xcos, suCOS, etc. Those scores are based on shape and chemical complementarity.
    """

    def __init__(self, fragments_dir=None, fragment_id_pattern=None, fragments_dict=None, selected_fragment_ids=None, *args, **kwargs):
        '''
        This params are generally provided through the class method computeScoreForMolecules directly obtained from cmd parser
        '''

        if not fragments_dict:
            assert fragments_dir is not None and fragment_id_pattern is not None, "Error, if no list of fragments (Chem.Mol) provided, " \
                                                                                  "you must specify the fragments_dir and fragment_id_pattern"
            fragments = load_files_as_mols( fragments_dir, file_pattern=fragment_id_pattern)
            fragments_dict = dict( fragments )
        else:
            assert fragments_dir is None and fragment_id_pattern is None, "Error, if a list of fragments (Chem.Mol) provided, " \
                                                                          "you should not specify fragments_dir and fragment_id_pattern"
        self.fragments_dict = fragments_dict
        if selected_fragment_ids:
            self.fragments_dict = { key: val for key,val in self.fragments_dict.items() if key in selected_fragment_ids}

        super().__init__( *args, **kwargs)

    @property
    def fragments_id(self):
        '''
        This is instantiation of abstract attribute
        :return:
        '''
        return list( list(self.fragments_dict.keys() ) )

    @classmethod
    def save_grid_as_mrc(cls, grid, outname, spacing=None):
        assert outname.endswith(outname), "Error, outname must be .mrc file"

        if isinstance(grid, rdGeometry.UniformGrid3D_):
            spacing = grid.GetSpacing() if spacing is None else spacing
            grid = cls.uniformGrid3D_to_numpy(grid)
        else:
            assert spacing is not None, "Error, if a numpy array used as argument, spacing is required"

        with mrcfile.new(outname, overwrite=True) as f:
            f.voxel_size = spacing
            f.set_data(grid.astype(np.float32))


def test1():

    #TODO: work only with the binding site, not the whole protein
    apo_pdb_mol = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x10889_0A/Mpro-x10889_0A_apo-desolv.pdb") )

    mol_0  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x10889_0A/Mpro-x10889_0A.mol") )
    mol_1  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x0755_0A/Mpro-x0755_0A.mol") )
    mol_2  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x0769_0A/Mpro-x0769_0A.mol") )
    mol_3  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x1308_0A/Mpro-x1308_0A.mol") )
    mol_4  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x1334_0A/Mpro-x1334_0A.mol") )
    mol_5  = load_mol( os.path.expanduser("./test_mols/Mpro-x1334_0A_translated.mol"))

    mols_list = [mol_0, mol_1, mol_2, mol_3, mol_4 ] + [mol_5]

    results_weight = _COSLikeBase.compute_occupancy_weights( mols_list )
    print( results_weight )
    results_clash = _COSLikeBase.compute_clashes( apo_pdb_mol, mols_list)
    print( results_clash )


def test2():

    mol_0 = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0176_0B/nsp13-x0176_0B.mol") )
    mol_1 = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0246_0B/nsp13-x0246_0B.mol") )
    mol_2 = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0176_0B/nsp13-x0176_0B.mol") )

    mols_list = [mol_0 , mol_1, mol_2 ]

    results_weight = _COSLikeBase.compute_occupancy_weights( mols_list )
    print( results_weight )
    results_weight = _COSLikeBase.compute_occupancy_weights( mols_list)
    print( results_weight )
    print( _COSLikeBase.computeSuCOS(mol_0, mol_1) )
    #TODO: check if makes sense

if __name__ == "__main__":
    test1()
    test2()

    '''

python -m fragmenstein.scoring.cos_like_base

    '''