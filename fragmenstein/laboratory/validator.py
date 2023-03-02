"""
This is an optional import and can only be accessed via `fragmenstein.laboratory.validator`.

For place a route is providing a dataframe with the following columns:

- ``name``: the name of the built molecule (_e.g._'build#👾')
- ``smiles``: the smiles of the built molecule
- ``hits``: the hits used to build the molecule

These can be passed to ``Laboratory.place`` to place the built molecules.
This is a pandera schema for validation.

.. code-block:: python
    from fragmenstein.laboratory.validator import place_input_validator
    place_input_validator.validate(df)

    from fragmenstein import Laboratory
    Laboratory(pdbblock=👾👾👾, covalent_resi=None).place(df, expand_isomers=False, n_cores=12)
"""

import pandas as pd
import pandera as pa
from rdkit import Chem


def hits_check(mols_col: pd.Series) -> bool:
    for mols in mols_col:
        if not isinstance(mols, list):
            return False
        for mol in mols:
            if not isinstance(mol, Chem.Mol):
                return False
    return True


def smiles_check(smiles_col: pd.Series) -> bool:
    for smiles in smiles_col:
        if Chem.MolFromSmiles(smiles) is None:
            return False
    return True


place_input_validator = pa.DataFrameSchema({
    "name": pa.Column(pa.String),
    "hits": pa.Column(checks=[pa.Check(hits_check)]),
    "smiles": pa.Column(pa.String, checks=[pa.Check(smiles_check)])})
