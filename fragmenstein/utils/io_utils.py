import os
import re
import dask #To provide parallel suport

from typing import Union
from rdkit import Chem

from fragmenstein.utils.config_manager import ConfigManager


def load_mol(fname: str) -> Chem.Mol:
    '''
    Loads a pdb or mol file as a rdkit.Chem.Mol

    :param fname: a pdb fname or mol fname
    :return: rdkit.Chem.Mol
    '''
    if fname.endswith(".pdb"):
        mol = Chem.MolFromPDBFile( fname )
    elif fname.endswith(".mol"):
        mol = Chem.MolFromMolFile( fname )
    else:
        raise NotImplementedError( "Only .pdb and .mol formats are supported when loading %s" %fname)

    Chem.SanitizeMol(mol)
    return mol

#TODO: check code duplication with Victor and moster utils

def load_mol_if_str( mol_or_str: Union[str, Chem.Mol]) -> Union[Chem.Mol, None]:
    '''
    if input is Chem.Mol object, do nothing, elseif string, try to load as Chem.Mol
    :param mol_or_str:
    :return: The same Chem.Mol as input if it was a Chem.Mol, or the Chem.Mol associated to input string or None
             if loading was unsuccessful.
    '''
    if isinstance(mol_or_str, str):
        try:
            return load_mol(mol_or_str)
        except OSError:
            print("OSError", mol_or_str, )
            return None
    elif isinstance(mol_or_str, Chem.Mol):
        return mol_or_str
    else:
        raise ValueError("Expected Chem.Mol or str. Found: %s "%type(mol_or_str))


def apply_func_to_files( folder, file_pattern, function, use_parallel_dask=None, extensions_to_check=None):
        results = []

        if use_parallel_dask is None:
            use_parallel_dask = ConfigManager.N_CPUS > 0

        if use_parallel_dask:
            function = dask.delayed(function)
        for root, dirs, fnames in os.walk(folder):
            for fname in fnames:
                match_obj = re.match(file_pattern, fname)
                if match_obj:
                    if extensions_to_check and not os.path.splitext(fname)[-1] in extensions_to_check:
                        raise ValueError(
                            "Error, one of the files found (%s) did not match the required extension (%s)" % (
                            fname, extensions_to_check))
                    fname = os.path.join( root, fname)
                    result = function(fname)
                    results.append(result)

        if use_parallel_dask:
            return  dask.compute(results)[0]
        else:
            return results


def load_files_as_mols( mols_folder, file_pattern=".*-(\w+_\w{2}).mol$", use_parallel_dask=None,
                       extensions_to_check=None):
    def process(fname):
        mol_id = re.match(file_pattern, os.path.split(fname)[-1]).group(1)
        mol = load_mol(fname)
        return mol_id, mol

    return apply_func_to_files(mols_folder, file_pattern=file_pattern, function=process,
                               use_parallel_dask=use_parallel_dask, extensions_to_check=extensions_to_check)
