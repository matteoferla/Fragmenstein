#!/usr/bin/env python

"""
This is a base clase for the COS-like scores such as xcos, suCOS, etc. Those scores are based on shape and chemical complementarity.
"""

import os
from copy import deepcopy
from functools import reduce

import dask
import joblib
import mrcfile
from dask.distributed import progress
import dask.bag as db
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit import RDConfig

import numpy as np
from rdkit.Chem.rdMolTransforms import ComputeCentroid, TransformConformer
from rdkit.Chem.rdShapeHelpers import ComputeUnionBox, ComputeConfDimsAndOffset, ComputeConfBox, EncodeShape
from rdkit.Geometry import rdGeometry

from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.parallel_utils import get_parallel_client

from fragmenstein.scoring._scorer_base import _ScorerBase, journal
from fragmenstein.utils.io_utils import load_mol, load_files_as_mols


def uniformGrid3D_to_numpy( grid):

    numpy_result = np.zeros( grid.GetSize())
    DataStructs.ConvertToNumpyArray(grid.GetOccupancyVect(), numpy_result)
    numpy_result = numpy_result.reshape((grid.GetNumZ(), grid.GetNumY(), grid.GetNumX()))
    return numpy_result

def compute_mol_to_grid(mol, grid_config, as_numpy):
    grid = rdGeometry.UniformGrid3D(*grid_config["size"], spacing=grid_config["spacing"])
    mol = deepcopy(mol)
    tfm= np.eye(4)
    tfm[:3,-1] =  - grid_config["coordinates_origin"]
    AllChem.TransformMol(mol, tfm)
    Chem.rdShapeHelpers.EncodeShape(mol, grid, ignoreHs=False, vdwScale=grid_config["vdwScale"])
    # save_grid_as_mrc(grid, "prueba_%d.mrc"%i)

    if as_numpy:
        grid = uniformGrid3D_to_numpy(grid)
    return grid


def save_grid_as_mrc( grid, outname, spacing=None):
    assert outname.endswith( outname), "Error, outname must be .mrc file"

    if isinstance( grid, rdGeometry.UniformGrid3D_ ):
        spacing = grid.GetSpacing() if spacing is None else spacing
        grid = uniformGrid3D_to_numpy(grid)
    else:
        assert  spacing is not None, "Error, if a numpy array used as argument, spacing is required"

    with mrcfile.new(outname,overwrite=True) as f:
        f.voxel_size = spacing
        f.set_data( grid.astype(np.float32))


_FEATURES_FACTORY=[]
def get_features_factory():
    if len(_FEATURES_FACTORY) ==0:
        # print("\n feature factory \n")
        feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
        fmParams = {k: FeatMaps.FeatMapParams() for k in feature_factory.GetFeatureFamilies()}
        keep_featnames = list(fmParams.keys())
        _FEATURES_FACTORY.extend([feature_factory, fmParams, keep_featnames])
    return _FEATURES_FACTORY

class _COSLikeBase(_ScorerBase):
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


    def getFeatureMapScore(self, small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.All):


        try:
            # feature_factory, fmParams, keep_featnames = self.feature_factory
            feature_factory, fmParams, keep_featnames = get_features_factory()
        except ValueError as e: #TODO: check that lock prevents race conditions and remove try/except
            raise e
        try:
            featLists = []
            for m in [small_m, large_m]:
                rawFeats = feature_factory.GetFeaturesForMol(m)
                # filter that list down to only include the ones we're intereted in
                featLists.append([f for f in rawFeats if f.GetFamily() in keep_featnames])
            fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params= fmParams) for x in featLists]
            fms[0].scoreMode = score_mode
            fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
            return fm_score

        except RuntimeError:
            return np.nan

        except ZeroDivisionError:
            return 0


    def getShapeScore(self, small_m, large_m, vdwScale=0.8):
        try:
            protrude_dist = rdShapeHelpers.ShapeProtrudeDist(small_m, large_m, allowReordering=False, vdwScale=vdwScale)
        except RuntimeError:
            return np.nan

        protrude_dist = np.clip(protrude_dist, 0, 1)
        protrude_score = 1 - protrude_dist
        return protrude_score

    @classmethod
    def estimate_grid_params(cls, list_of_mols, spacing=0.4, vdwScale=0.8, box_margin=4, expansion_factor=2):

        get_BoxLimits = lambda mol : Chem.rdShapeHelpers.ComputeConfBox( mol.GetConformer())
        box_limits = get_BoxLimits(list_of_mols[0])
        coordinates_origin = ComputeCentroid(list_of_mols[0].GetConformer())

        for mol in list_of_mols[1:]:
            box_limits = ComputeUnionBox(box_limits, get_BoxLimits(mol))

        size = (box_limits[1].x - box_limits[0].x + 2 * box_margin,
                box_limits[1].y - box_limits[0].y + 2 * box_margin,
                box_limits[1].z - box_limits[0].z + 2 * box_margin)
        size = tuple(expansion_factor * num for num in size)

        grid_config = {"size":size, "spacing":spacing, "vdwScale":vdwScale, "box_margin":box_margin,
                       "expansion_factor":expansion_factor, "coordinates_origin":np.array([coordinates_origin.x,
                                                                           coordinates_origin.y, coordinates_origin.z])}
        return grid_config

    @classmethod
    def yield_uniformGrid3D(cls, list_of_mols, grid_config=None, spacing=0.4, vdwScale=0.8, box_margin=4, expansion_factor=2,
                            as_numpy=False):
        if grid_config is None:
            grid_config = cls.estimate_grid_params( list_of_mols, spacing, vdwScale, box_margin, expansion_factor)
        for mol in list_of_mols:
            grid = compute_mol_to_grid(mol, grid_config, as_numpy)
            yield grid

    @classmethod
    def compute_occupancy_weights(cls, mols_list, vdwScale=0.2, spacing=0.3, use_joblib_instead_dask=False):

        # print(len(mols_list)); mols_list = mols_list[:30]; print("WARNING, debug MODE")

        journal.info( "Loading molecules as grids to compute weights")
        print( "Loading molecules as grids to compute weights")

        grid_config = cls.estimate_grid_params(mols_list, spacing, vdwScale)

        def process_one_grid(mol, sumGrid, sumGrid_nonZero):
            grid = compute_mol_to_grid(mol, grid_config, as_numpy=True)
            n_vox_mol =  np.count_nonzero(grid)
            grid[sumGrid_nonZero] /= sumGrid[sumGrid_nonZero]
            result =  np.sum(grid) / n_vox_mol
            return result

        if use_joblib_instead_dask:
            all_grids = joblib.Parallel(n_jobs = ConfigManager.N_CPUS)( joblib.delayed(compute_mol_to_grid)(mol, grid_config, True) for mol in mols_list )
            sumGrid = reduce(lambda prev, cur : prev + cur, all_grids, 0)
            sumGrid_nonZero = ~ np.isclose(sumGrid, 0)
            final_weights = [process_one_grid(mol, sumGrid, sumGrid_nonZero) for mol in  mols_list ]

        else:
            client = get_parallel_client()
            b = db.from_sequence( mols_list )
            mapped_b = b.map(lambda mol: compute_mol_to_grid(mol, grid_config, True)).fold(lambda prev, cur : prev + cur, initial=0)
            sumGrid_future = client.compute(mapped_b)

            sumGrid_future = client.compute(sumGrid_future)
            sumGrid_nonZero_future = client.submit( lambda x: ~ np.isclose(x, 0), sumGrid_future)
            del b
            final_weights = map( lambda mol: client.submit(process_one_grid, *(mol, sumGrid_future, sumGrid_nonZero_future)), mols_list)
            final_weights = client.compute(final_weights)
            final_weights= client.gather(final_weights)

        journal.info( "Grid weights computed")
        # print(final_weights)
        return final_weights


    @classmethod
    def compute_clashes(cls, ref_mol, mols_list, vdwScale=0.4, spacing=0.5):
        '''

        :param ref_mol:
        :param mol_list:
        :param vdwScale:
        :param spacing:
        :return:
        '''


        mols_list = list(mols_list)
        grid_config = _COSLikeBase.estimate_grid_params([ref_mol] + mols_list, vdwScale=vdwScale, spacing=spacing)

        if isinstance(ref_mol, Chem.Mol):
            ref_grid = compute_mol_to_grid(ref_mol, grid_config=grid_config, as_numpy=True)
        else:
            ref_grid = ref_mol

        if isinstance(mols_list[0], Chem.Mol):
            mols_grids = [dask.delayed(compute_mol_to_grid)(mol, grid_config, as_numpy=True) for mol in mols_list]
        else:
            mols_grids = mols_list

        def computeClash(grid):
            grid_intersect = (grid * ref_grid)
            return np.sum(grid_intersect) / np.sum(grid)  # np.count_nonzero

        clashes_list = []
        for i, grid in enumerate(mols_grids):
            score = dask.delayed(computeClash)(grid)
            clashes_list.append(score)

        client = get_parallel_client()
        clashes_list = client.compute(clashes_list)
        progress( clashes_list )
        clashes_list = client.gather(clashes_list)

        return clashes_list


    # @classmethod
    # def compute_clash(cls, ref_mol, mols_list, vdwScale=0.4, spacing=0.5):
    #     '''
    #
    #     :param ref_mol:
    #     :param mol_list:
    #     :param vdwScale:
    #     :param spacing:
    #     :return:
    #     '''
    #     mols_list = list( mols_list )
    #     grid_config = _COSLikeBase.estimate_grid_params([ref_mol] + mols_list, vdwScale=vdwScale, spacing=spacing)
    #
    #     if isinstance(ref_mol, Chem.Mol ):
    #         ref_grid = next( _COSLikeBase.yield_uniformGrid3D([ref_mol], grid_config=grid_config, as_numpy=True) )
    #     else:
    #         ref_grid = ref_mol
    #
    #     if isinstance(mols_list[0], Chem.Mol ):
    #         mols_grids = next( _COSLikeBase.yield_uniformGrid3D( mols_list, grid_config=grid_config, as_numpy=True) )
    #     else:
    #         mols_grids = mols_list
    #
    #     sumGrid =  0
    #     for grid in mols_grids:
    #         sumGrid += grid
    #
    #     if isinstance(mols_list[0], Chem.Mol ):
    #         mols_grids =  _COSLikeBase.yield_uniformGrid3D( mols_list, grid_config=grid_config, as_numpy=True)
    #     else:
    #         mols_grids = mols_list
    #
    #     clashes_list = [0]* len(mols_list)
    #     for i, grid in enumerate(mols_grids):
    #         grid_intersect = (grid * ref_grid)
    #         clashes_list[i] = np.sum(grid_intersect) / np.sum( grid) # np.count_nonzero
    #
    #     return clashes_list


def test1():
    dask_client = get_parallel_client(); print(dask_client)

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
    results_weight = _COSLikeBase.compute_occupancy_weights( mols_list , use_joblib_instead_dask=True)
    print( results_weight )
    #TODO: check if makes sense

if __name__ == "__main__":
    # test1()
    test2()

    '''

python -m fragmenstein.scoring.cos_like_base

    '''