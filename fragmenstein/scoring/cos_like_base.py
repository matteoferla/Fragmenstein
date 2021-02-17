#!/usr/bin/env python

"""
This is a base clase for the COS-like scores such as xcos, suCOS, etc. Those scores are based on shape and chemical complementarity.
"""

import os
import threading
import dask
import mrcfile
from dask.distributed import progress
from functools import reduce
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit import RDConfig

import numpy as np
from rdkit.Chem.rdShapeHelpers import ComputeUnionBox
from rdkit.Geometry import rdGeometry
from fragmenstein.utils.parallel_utils import get_parallel_client

from fragmenstein.scoring._scorer_base import _ScorerBase, journal
from fragmenstein.utils.io_utils import apply_func_to_files, load_mol, load_files_as_mols

feature_factory_objs = [] #this global object used as a cache is required because BuildFeatureFactory produces non pickleable objects
def get_feature_factory():
    if len(feature_factory_objs)==0:
        feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
        fmParams = {k: FeatMaps.FeatMapParams() for k in feature_factory.GetFeatureFamilies()}
        keep_featnames = list(fmParams.keys())
        feature_factory_objs.extend([feature_factory, fmParams, keep_featnames ] )
    return feature_factory_objs[:3]




def uniformGrid3D_to_numpy( grid):

    numpy_result = np.zeros( grid.GetSize())
    DataStructs.ConvertToNumpyArray(grid.GetOccupancyVect(), numpy_result)
    numpy_result = numpy_result.reshape((grid.GetNumZ(), grid.GetNumY(), grid.GetNumX()))
    return numpy_result

def compute_one_grid( mol, grid_config, as_numpy):
    grid = rdGeometry.UniformGrid3D(*grid_config["size"], spacing=grid_config["spacing"])
    Chem.rdShapeHelpers.EncodeShape(mol, grid, ignoreHs=False, vdwScale=grid_config["vdwScale"])
    if as_numpy:
        grid = uniformGrid3D_to_numpy(grid)
    return grid


class _COSLikeBase(_ScorerBase):
    """
    This is a base clase for the COS-like scores such as xcos, suCOS, etc. Those scores are based on shape and chemical complementarity.
    """

    def __init__(self, fragments_dir, fragment_id_pattern, *args, **kwargs):
        '''
        This params are generally provided through the class method computeScoreForMolecules directly obtained from cmd parser
        '''


        fragments = load_files_as_mols( fragments_dir, file_pattern=fragment_id_pattern)
        self.fragments_dict = dict(fragments)

        super().__init__( *args, **kwargs)

        try:
            self.lock = dask.distributed.Lock("feature_factory_objs")
        except ValueError:
            self.lock = threading.Lock()

    @property
    def fragments_id(self):
        '''
        This is instantiation of abstract attribute
        :return:
        '''
        return list( list(self.fragments_dict.keys() ) )


    def getFeatureMapScore(self, small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.All):


        try:
            with self.lock:
                feature_factory, fmParams, keep_featnames = get_feature_factory()

        except ValueError as e: #TODO: check that lock prevents race conditions and remove try/except
            print(get_feature_factory())
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

        for mol in list_of_mols[1:]:
            box_limits = ComputeUnionBox(box_limits, get_BoxLimits(mol))

        size = (box_limits[1].x - box_limits[0].x + 2 * box_margin,
                box_limits[1].y - box_limits[0].y + 2 * box_margin,
                box_limits[1].z - box_limits[0].z + 2 * box_margin)
        size = tuple(expansion_factor * num for num in size)

        grid_config = {"size":size, "spacing":spacing, "vdwScale":vdwScale, "box_margin":box_margin,
                       "expansion_factor":expansion_factor}
        return grid_config

    @classmethod
    def yield_uniformGrid3D(cls, list_of_mols, grid_config=None, spacing=0.4, vdwScale=0.8, box_margin=4, expansion_factor=2,
                            as_numpy=False):
        if grid_config is None:
            grid_config = cls.estimate_grid_params( list_of_mols, spacing, vdwScale, box_margin, expansion_factor)
        for mol in list_of_mols:
            grid = compute_one_grid(mol, grid_config, as_numpy)
            yield grid


    @classmethod
    def save_grid_as_mrc(cls, grid, outname, spacing=None):

        assert outname.endswith( outname), "Error, outname must be .mrc file"

        if isinstance( grid, rdGeometry.UniformGrid3D_ ):
            spacing = grid.GetSpacing() if spacing is None else spacing
            grid = uniformGrid3D_to_numpy(grid)
        else:
            assert  spacing is not None, "Error, if a numpy array used as argument, spacing is required"

        with mrcfile.new(outname,overwrite=True) as f:
            f.voxel_size = spacing
            f.set_data( grid.astype(np.float32))

    @classmethod
    def compute_occupancy_weights(cls, mols_list, vdwScale=0.2, spacing=0.3):

        # print(len(mols_list)); mols_list = mols_list[:30]; print("WARNING, debug MODE")

        journal.info( "Loading molecules as grids to compute weights")
        print( "Loading molecules as grids to compute weights")

        grid_config = cls.estimate_grid_params(mols_list, spacing, vdwScale)

        client = get_parallel_client()

        #This option seems to eat much more memory
        # sumGrid = 0
        # for mol in mols_list:
        #     grid =  dask.delayed(compute_one_grid)( mol, grid_config, as_numpy=True)
        #     sumGrid += grid
        # sumGrid = client.compute(sumGrid)
        # progress( sumGrid )
        # sumGrid= sumGrid.result()

        grids_delayed_gen = ( dask.delayed(compute_one_grid)( mol, grid_config, True) for mol in mols_list)
        sumGrid_future = reduce(lambda prev, cur : prev+cur, grids_delayed_gen)
        # as_completed()

        sumGrid_future = client.compute(sumGrid_future)
        sumGrid_nonZero_future = client.submit( lambda x: ~ np.isclose(x, 0), sumGrid_future)

        del grids_delayed_gen

        def process_one_grid(mol, sumGrid, sumGrid_nonZero):
            grid = compute_one_grid( mol, grid_config, as_numpy=True)
            n_vox_mol =  np.count_nonzero(grid)
            grid[sumGrid_nonZero] /= sumGrid[sumGrid_nonZero]
            result =  np.sum(grid) / n_vox_mol
            return result

        final_weights = map( lambda mol: client.submit(process_one_grid, *(mol, sumGrid_future, sumGrid_nonZero_future)), mols_list)

        final_weights = client.compute(final_weights)
        progress( final_weights )
        final_weights= client.gather(final_weights)

        journal.info( "\nWeights computed")
        print( "\nWeights computed")
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
            ref_grid = _COSLikeBase.compute_one_grid(ref_mol, grid_config=grid_config, as_numpy=True)
        else:
            ref_grid = ref_mol

        if isinstance(mols_list[0], Chem.Mol):
            mols_grids = [ dask.delayed(compute_one_grid)( mol, grid_config, as_numpy=True) for mol in  mols_list]

            def computeClash(mol):
                grid = compute_one_grid( mol, grid_config, as_numpy=True)
                grid_intersect = (grid * ref_grid)
                return np.sum(grid_intersect) / np.sum(grid)  # np.count_nonzero

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


    @classmethod
    def compute_clash(cls, ref_mol, mols_list, vdwScale=0.4, spacing=0.5):
        '''

        :param ref_mol:
        :param mol_list:
        :param vdwScale:
        :param spacing:
        :return:
        '''
        mols_list = list( mols_list )
        grid_config = _COSLikeBase.estimate_grid_params([ref_mol] + mols_list, vdwScale=vdwScale, spacing=spacing)

        if isinstance(ref_mol, Chem.Mol ):
            ref_grid = next( _COSLikeBase.yield_uniformGrid3D([ref_mol], grid_config=grid_config, as_numpy=True) )
        else:
            ref_grid = ref_mol

        if isinstance(mols_list[0], Chem.Mol ):
            mols_grids = next( _COSLikeBase.yield_uniformGrid3D( mols_list, grid_config=grid_config, as_numpy=True) )
        else:
            mols_grids = mols_list

        sumGrid =  0
        for grid in mols_grids:
            sumGrid += grid

        if isinstance(mols_list[0], Chem.Mol ):
            mols_grids =  _COSLikeBase.yield_uniformGrid3D( mols_list, grid_config=grid_config, as_numpy=True)
        else:
            mols_grids = mols_list

        clashes_list = [0]* len(mols_list)
        for i, grid in enumerate(mols_grids):
            grid_intersect = (grid * ref_grid)
            clashes_list[i] = np.sum(grid_intersect) / np.sum( grid) # np.count_nonzero

        return clashes_list


if __name__ == "__main__":
    dask_client = get_parallel_client(); print(dask_client)

    #TODO: work only with the binding site, not the whole protein
    apo_pdb_mol = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x10889_0A/Mpro-x10889_0A_apo-desolv.pdb") )
    mol_0  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x10889_0A/Mpro-x10889_0A.mol") )
    mol_1  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x0755_0A/Mpro-x0755_0A.mol") )
    mol_2  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x0769_0A/Mpro-x0769_0A.mol") )
    mol_3  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x1308_0A/Mpro-x1308_0A.mol") )
    mol_4  = load_mol( os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro/aligned/Mpro-x1334_0A/Mpro-x1334_0A.mol") )
    mol_5  = load_mol( os.path.expanduser("~/tmp/Mpro-x1334_0A_translated.mol"))

    mols_list = [mol_0, mol_1, mol_2, mol_3, mol_4 ] + [mol_5]

    results_weight = _COSLikeBase.compute_occupancy_weights( mols_list )
    print( results_weight )
    results_clash = _COSLikeBase.compute_clash( apo_pdb_mol, mols_list)
    print( results_clash )

    '''

python -m fragmenstein.scoring.cos_like_base

    '''