from rdkit import Chem
from rdkit.Chem import BRICS, AllChem
from typing import List, Tuple, Any, Dict
import numpy as np
import pandas as pd

class AccountableBRICS:
    """
    BRICS decomposition that remembers where stuff came from.

    This stores the original inspiration as an isotope label
    It is potentially dangerous as dummy atom isotope number plays a role too.
    No property is kept however, so it is the only working solution.
    The attachment atom is labelled in the built molecule with `bridging_atom` which contains its inspiration's name

    ... code-block:: python
        decomposer = AccountableBRICS(hits)
        df: pd.DataFrame = decomposer(cutoff=decomposer.median, max_mergers=100)
        decomposer.info
        # {'N_hits': 44, 'N_fragments': 157, 'too_small': 0, 'missing': [], 'max_reached': True, 'N_built': 100}

    The output of ``__call__`` is a pandas dataframe with the following columns:

    - ``name``: the name of the built molecule ('build#👾')
    - ``smiles``: the smiles of the built molecule
    - ``hits``: the hits used to build the molecule
    - ``built_molecule``: the built molecule

    These can be passed to ``Laboratory.place`` to place the built molecules.

    ... code-block:: python
        pyrosetta.init(extra_options=extra_options)
        lab = Laboratory(pdbblock=👾👾👾, covalent_resi=None)
        placements:pd.DataFrame = lab.place(df, expand_isomers=False, n_cores=12)
        display(placements)
    """

    def __init__(self, hits: List[Chem.Mol]):
        self.hits = {}
        for i, hit in enumerate(hits):
            hit.SetIntProp('decompose_index', i + 1000)
            self.hits[i + 1000] = hit
        self.info = {'N_hits': len(self.hits),
                     'N_fragments': 0,
                     'too_small': 0,
                     'missing': [],
                     'max_reached': False}
        self.fragments: List[Chem.Mol] = []

    def decompose(self, mol: Chem.Mol) -> List[Chem.Mol]:
        name = mol.GetProp('_Name')
        index = mol.GetIntProp('decompose_index')
        for atom in mol.GetAtoms():
            atom.SetProp('ori_name', name)
            atom.SetIntProp('ori_i', atom.GetIdx())
        fragments = list(BRICS.BRICSDecompose(mol, keepNonLeafNodes=False, returnMols=True))
        # label
        for i, frag in enumerate(fragments):
            # this is a waste of time:
            frag.SetProp('_Name', f'{name}#{i}')
            for atom in frag.GetAtoms():
                atom.SetProp('ori_name', name)
            # but isotope is kept:
            dummy: Chem.Atom
            for dummy in frag.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(0)):
                dummy.GetNeighbors()[0].SetIsotope(index)
                pass
        return fragments

    def __call__(self, cutoff: int = 0, max_mergers: int = 10_000, maxDepth: int = 4) -> pd.DataFrame:
        self.fragments = []
        for hit in self.hits.values():
            self.fragments.extend(self.decompose(hit))
        self.info['N_fragments'] = len(self.fragments)
        mergers: List[Chem.Mol] = []
        accepted: int = 0
        results: List[Dict[str, Any]] = []
        for built in BRICS.BRICSBuild(self.fragments,
                                      onlyCompleteMols=True,
                                      maxDepth=maxDepth):
            if built.GetNumHeavyAtoms() < cutoff:
                self.info['too_small'] += 1
                continue
            accepted += 1
            inspirations: List[Chem.Mol] = self.get_inspirations(built)
            results.append({'name': f'build#{len(results)}',  # required by Laboratory
                            'built_molecule': built,
                            'smiles': Chem.MolToSmiles(built),  # required by Laboratory
                            'hits': inspirations,
                            })
            if accepted >= max_mergers:
                self.info['max_reached'] = True
                break
        self.info['N_built'] = len(results)
        return pd.DataFrame(results)

    @property
    def median(self):
        return np.median(list(map(Chem.Mol.GetNumHeavyAtoms, self.hits.values())))

    def get_inspirations(self, built: Chem.Mol) -> List[Chem.Mol]:
        inspirations: List[Chem.Mol] = []
        atom: Chem.Atom
        for atom in built.GetAtomsMatchingQuery(AllChem.IsotopeGreaterQueryAtom(0)):
            i = atom.GetIsotope()
            if i not in self.hits:
                self.info['missing'].append(i)
                continue
            inspirations.append(self.hits[i])
            atom.SetIsotope(0)
            atom.SetProp('bridging_atom', self.hits[i].GetProp('_Name'))
        return sorted(inspirations, key=Chem.Mol.GetNumHeavyAtoms, reverse=True)