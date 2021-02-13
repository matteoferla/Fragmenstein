## INCOMPLETE

raise NotImplementedError


from sqlitedict import SqliteDict
from rdkit.Chem import PandasTools
import json
import pandas as pd
from fragmenstein.victor import Victor

import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import skewnorm, gennorm
from typing import Dict, List, Union

from ._process import process

class Laboratory:

    def __init__(self, project, hits: list[Chem.Mol]):
        self.project = project
        self.hits = {hit.GetProp('_Name'): hit for hit in hits}
        if None in self.hits:
            raise ValueError('Molecule without a name given.')

    def merge(self, cores = 25):
        ## Process
        pass






    def old_ranker(self, row):
        try:
            return float(row['∆∆G']) / 5 + float(
                row.comRMSD) + row.N_unconstrained_atoms / 5 - row.N_constrained_atoms / 10
            # return float(row['∆∆G'])/(row.N_unconstrained_atoms + row.N_constrained_atoms * 0.5)*10 + float(row.comRMSD)
        except:
            return float('nan')

    rank_weights = {'LE': 1., 'comRMSD': 2., 'atom_bonus': 2., 'novelty_penalty': 5.}

    def ranker(self, row):
        try:
            # atom_bonus = row.N_constrained_atoms / (20 + row.N_constrained_atoms)
            # atom_bonus = skewnorm.pdf((row.N_constrained_atoms - 20)/8, 3)
            ζ = (row.N_constrained_atoms ** 2 - 25 ** 2) / 500
            atom_bonus = gennorm.pdf(ζ, 5) / 0.5445622105291682
            novelty_penalty = row.N_unconstrained_atoms / row.N_constrained_atoms
            return rank_weights['LE'] * float(row.LE) + \
                   rank_weights['comRMSD'] * float(row.comRMSD) + \
                   - rank_weights['atom_bonus'] * atom_bonus + \
                   rank_weights['novelty_penalty'] * novelty_penalty
        except:
            return float('nan')

    def LE(self, row):
        try:
            return float(row['∆∆G']) / (row.N_unconstrained_atoms + row.N_constrained_atoms)
        except:
            return float('nan')

    def get_mol3D(self, name):
        path = os.path.join(Victor.work_path, name, name + '.minimised.mol')
        if os.path.exists(path):
            try:
                mol = Chem.MolFromMolFile(path, sanitize=True)
                if mol is None:
                    return None
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                return mol
            except Exception as error:
                print(f'{type(error)}: {error}')
                return None
        else:
            return None

    def get_table(self, db_name, mols=True, mol_only=True):
        results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
        result_table = pd.DataFrame(results.values())
        print(len(result_table), sum(~result_table['∆∆G'].isna()))
        result_table['LE'] = result_table.apply(LE, 1)
        rank = result_table.apply(ranker, axis=1).rank()
        m = np.nanmax(rank.values)
        result_table['%Rank'] = rank / m * 100
        result_table['N_hits'] = result_table.regarded.apply(lambda x: len(x) if str(x) != 'nan' else float('nan'))
        result_table = result_table.loc[~result_table.smiles.isna()].sort_values(['%Rank'], axis=0)
        if mols:
            result_table['mol3D'] = result_table['name'].apply(get_mol3D)
            # result_table['mol2D'] = result_table['name'].apply(get_mol2D)
            PandasTools.AddMoleculeColumnToFrame(result_table, 'smiles', 'mol2D')
            if mol_only:
                result_table = result_table.loc[~result_table.mol3D.isna()]
        return result_table

    atom_Ns = {}

    for folder in ('newinputs',):  # 'input', 'UCSF2-hits', 'frags'):
        for file in os.listdir(folder):
            if '.mol' in file:
                mol = Chem.MolFromMolFile(os.path.join(folder, file), sanitize=False)
                if mol is None:
                    atom_Ns[file.replace('.mol', '')] = 0  # float nan?
                else:
                    mol = Chem.GetMolFrags(mol, asMols=True)[0]  # just in case
                    atom_Ns[file.replace('.mol', '')] = mol.GetNumAtoms()