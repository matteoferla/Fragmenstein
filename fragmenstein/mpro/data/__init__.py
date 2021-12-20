from . import hit_mols

import os, pyrosetta
from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import List, Optional

import pandas as pd
import io
import requests
import random
import importlib.resources as pkg_resources

def get_template():
    return pkg_resources.read_text(__package__, 'template.pdb')

def _clean_hitname(hit_name: Optional[str]=None) -> str:
    if hit_name is None:
        hit_name = random.choice(get_hit_list())
        'Mpro-6y2g.mol'
    hit_name = hit_name.replace('.mol', '').replace('Mpro-', '').strip()
    if not pkg_resources.is_resource(hit_mols, f'Mpro-{hit_name}.mol'):
        raise FileNotFoundError(f'There is no hit {hit_name} in the cache. Choices are {get_hit_list()}')
    return hit_name

def get_molblock(hit_name: Optional[str]=None) -> str:
    """
    returns the mol block for the given hit. Chosen randomly if unspecified.
    See ``get_mpro_hit_list()`` for options.
    """
    hit_name = _clean_hitname(hit_name)
    return pkg_resources.read_text(hit_mols, f'Mpro-{hit_name}.mol')

def get_mol(hit_name: Optional[str]=None) -> Chem.Mol:
    """
    returns the Chem.Mol instance for the given hit. Chosen randomly if unspecified.
    See ``get_mpro_hit_list()`` for options.
    """
    hit_name = _clean_hitname(hit_name)
    mol = Chem.MolFromMolBlock(get_molblock(hit_name))
    mol.SetProp('_Name', hit_name)
    return mol

def get_filtered_mol(mw_cutoff: float=float('inf')) -> Chem.Mol:
    # For now only MolWt
    choices = get_hit_list()
    random.shuffle(choices)
    for hit_name in choices:
        mol = get_mol(hit_name)
        if mw_cutoff >= Descriptors.MolWt(mol):
            return mol
    else:
        raise ValueError('No hit matches the specified values')

def get_n_filtered_mols(amount: int, **cutoffs) -> List[Chem.Mol]:
    mols = []
    mol_names = set()
    assert amount > 0, 'A zero amount does nothing.'
    for i in range(len(get_hit_list())):
        mol = get_filtered_mol(**cutoffs)
        mol_name = mol.GetProp('_Name')
        if mol_name in mol_names:
            continue
        mol_names.add(mol_name)
        mols.append(mol)
        if len(mols) == int(amount):
            return mols
    else:
        raise ValueError(f'There are only {len(mols)} molecules ({mol_names}) that match the criteria {cutoffs}')

def get_hit_list():
    """
    List of XChem hits of MPro from Fragalysis in Feb 2021.
    """
    return [fn.replace('.mol', '').replace('Mpro-', '') for fn in pkg_resources.contents(hit_mols) if '.mol' in fn]