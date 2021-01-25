import os
import re
import warnings
import dask
import numpy as np


from dask.distributed import Client
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem import BRICS

from fragmenstein.scoring.scoring_config import  N_THREADS_PER_WORKERS, N_WORKERS, TMP_DIR


DASK_CLIENT = None
DASK_CLIENT_READY =False
def prepare_paralell_execution(threads_per_worker=N_THREADS_PER_WORKERS, n_workers=N_WORKERS):
    global DASK_CLIENT, DASK_CLIENT_READY
    if DASK_CLIENT is None:
        dask.config.set({'temporary_directory': os.path.join(TMP_DIR, "dask")})
        DASK_CLIENT = Client(threads_per_worker=threads_per_worker, n_workers=n_workers)
        DASK_CLIENT_READY = True

    return DASK_CLIENT


def mol(args):
    pass


class _ScorerUtils():
    '''
    Class that provides some utilities
    '''

    LIGAND_TAG="LIG"
    SKIP_RESIDUES_LIST = ["DMS"]

    @classmethod
    def load_molecule(cls, fname):
        if fname.endswith(".pdb"):
            mol = Chem.MolFromPDBFile( fname )
        elif fname.endswith(".mol"):
            mol = Chem.MolFromMolFile( fname )
        else:
            raise NotImplementedError("Only .pdb and .mol format are supported when loading %s"%fname)

        Chem.SanitizeMol(mol)
        return mol


    @classmethod
    def split_mol_to_Bits(cls, mol):
        '''
        :param mol. Chem.Mol object to be broken up into fragments by breaking rotable bonds
        :return:  a list of Chem.Mol objects that reprsent the bits in which input mol was broken.
        '''
        # find the rotatable bonds
        bonds = mol.GetSubstructMatches(RotatableBondSmarts)

        bonds = [((x, y), (0, 0)) for x, y in bonds]
        p = BRICS.BreakBRICSBonds(mol, bonds=bonds)
        mols = [mol for mol in Chem.GetMolFrags(p, asMols=True)]
        return mols

    @classmethod
    def load_mol_if_str(cls, mol_or_str):
        if isinstance(mol_or_str, str):
            try:
                mol = cls.load_molecule( mol_or_str)
            except OSError:
                print("OSError", mol_or_str, )
                mol = None
        return mol

    @classmethod
    def clean_pdb_mol(cls, pdb_mol, ):
        bad_residues = Chem.rdmolops.SplitMolByPDBResidues(pdb_mol, whiteList=_ScorerUtils.SKIP_RESIDUES_LIST, negateList=False)
        for br in bad_residues.values():
            pdb_mol = Chem.rdmolops.DeleteSubstructs(pdb_mol, br)
        return pdb_mol

    @classmethod
    def split_pdbMol_to_prot_lig(cls, pdb_mol, lig_mol=None, mol_id= None): #TODO; This may not work for covalently bonded structures

        warnings.warn("split_pdbMol_to_prot_lig is a experimental feature that does not work well for some covalent ligands" )
        pdb_mol = cls.clean_pdb_mol(pdb_mol)
        if lig_mol is not None:
            protein = Chem.rdmolops.DeleteSubstructs(pdb_mol, lig_mol)
            ligand = lig_mol

            mw_ligand = ExactMolWt(ligand)
            mw_boundPdb = ExactMolWt(pdb_mol)
            mw_protein = ExactMolWt(protein)

            if np.isclose( mw_protein, mw_boundPdb, atol=2.):
                raise Exception("Substructure was not removed")

            if mol_id:
              print( mol_id, mw_ligand, mw_boundPdb, mw_protein)

        else:
            try:
                ligand = Chem.rdmolops.SplitMolByPDBResidues(pdb_mol, whiteList=[_ScorerUtils.LIGAND_TAG])[_ScorerUtils.LIGAND_TAG]
                protein = Chem.rdmolops.DeleteSubstructs(pdb_mol, ligand)
            except KeyError as e:
                raise e

        return protein, ligand


    @classmethod
    def apply_func_to_files(cls, folder, file_pattern, function, delay_loading=False):
        results = []
        for root, dirs, fnames in os.walk(folder):
            for fname in fnames:
                match_obj = re.match(file_pattern, fname)
                if match_obj:
                    fname = os.path.join( root, fname)
                    result = dask.delayed(function)(fname)
                    results.append(result)

        if delay_loading:
            return results
        else:
            return dask.compute(results)[0]


    @classmethod
    def load_files_as_mols(cls, mols_folder, file_pattern="Mpro-(\w+_\w{2})_bound.pdb$", delay_loading=False):

        def process(fname):
            mol_id = re.match(file_pattern, os.path.split(fname)[-1]).group(1)
            mol = _ScorerUtils.load_molecule(fname)
            # mol.SetProp("_Name", mol_id)
            return mol_id, mol

        return cls.apply_func_to_files(mols_folder, file_pattern=file_pattern, function=process, delay_loading=delay_loading)

    @classmethod
    def load_pdbAndligand_as_mols(cls, results_folder, shared_file_pattern="Mpro-(\w+_\w{2})"):

        # assert DASK_CLIENT_READY is not None, "Error, Dask cliente should be first initialized using prepare_paralell_execution()"

        idAndPdb_mols = cls.load_files_as_mols(results_folder, file_pattern=shared_file_pattern + "_apo-desolv.pdb$", delay_loading=False)
        idAndLigand_mols = cls.load_files_as_mols(results_folder, file_pattern=shared_file_pattern+".mol$", delay_loading=False)

        results= {}
        for (mol_id, pdb), (_,ligand) in zip(idAndPdb_mols, idAndLigand_mols):
            # results.append( [mol_id, [pdb, ligand] ])results
            results[mol_id] = [pdb, ligand]
        return results



if __name__ == "__main__":

    from fragmenstein.scoring.scoring_config import MPRO_RAW_DATA_DIR

    desk_client = prepare_paralell_execution()

    results = _ScorerUtils.load_pdbAndligand_as_mols(os.path.join(MPRO_RAW_DATA_DIR, "aligned"))

    for mol_id, (pdb, lig) in results:
        print (mol_id, pdb, lig )
        Chem.MolToPDBFile(pdb, os.path.join(TMP_DIR, "test_dask", mol_id+"_prot.pdb"))
        Chem.MolToPDBFile(lig, os.path.join(TMP_DIR, "test_dask", mol_id+ "_lig.pdb"), )

