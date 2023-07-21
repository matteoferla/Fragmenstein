## This is an optional import

"""
Copy of
https://github.com/matteoferla/PLIP-PyRosetta-hotspots-test/blob/main/plipspots_docking/plipspots/serial.py

This is a class I use to _apply_ PLIP to a pd.Series of molecules.
It is not built for the project, but works.
Note I have not dumped the methods that are not needed for the project.
"""

import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools

from functools import singledispatchmethod
from typing import Tuple, Dict, List, Union, Sequence, Optional, Any
from collections import Counter, defaultdict
from plip.structure.preparation import PDBComplex, PLInteraction
from openbabel.pybel import Atom, Residue
from openbabel.pybel import ob
from .minimalPDB import MinimalPDB
import warnings


class SerialPLIPper:
    """
    Calling the instance will return a ``Dict[Tuple[str, str, int], int]``,
    where the key is interaction type, residue 3-letter name, residue index
    and the value is the count of interactions.
    Basically, applying Plip to a pd.Series of Chem.Mol.

    Unplacking it is kind of wierd, the best way I reckon is a brutal for-loop:

    .. code-block:: python

        import pandas as pd
        import pandera.typing as pdt

        intxndexes: pdt.Series[Dict[Tuple[str, str, int], int]] = hits.ROMol.apply(SerialPLIPper(pdb_filename))
        # columns will still be a tuple...:
        intxn_df = pd.DataFrame(intxndexes.to_list()).fillna(0).astype(int)
        hits['N_interactions'] = intxn_df.sum(axis='columns')
        for c in sorted(intxn_df.columns, key=lambda kv: kv[2]):
            # columns will be a colon-separated string:
            hits[':'.join(map(str, c))] = intxn_df[c]
    """

    def __init__(self, pdb_block: str, resn='LIG', resi=1, chain='B'):
        assert 'ATOM' in pdb_block, f'No ATOM entry in block provided: {pdb_block}'
        self.pdb_block = pdb_block
        self.resn = resn
        self.resi = int(resi)  # dont give me alt codes
        self.chain = chain

    @classmethod
    def from_filename(cls, pdb_filename: str, *args, **kwargs):
        """
        The main constructor is from PDB block, this is from PDB file
        """
        with open(pdb_filename, 'r') as f:
            pdb_block = f.read()
        return cls(pdb_block, *args, **kwargs)

    def __call__(self, mol) -> Dict[Tuple[str, str, int], int]:
        if mol is None or not isinstance(mol, Chem.Mol) or mol.GetNumAtoms() == 0:
            return {}
        holo: str = self.plonk(mol)
        interaction_set: PLInteraction = self.get_interaction_set(holo)
        return self.get_interaction_counts(interaction_set)

    def assign_pdb(self, mol: Chem.Mol):
        """
        Fix the PDB info for the molecule, in place
        """
        counts = defaultdict(int)
        atom: Chem.Atom
        for atom in mol.GetAtoms():
            element: str = atom.GetSymbol()
            counts[element] += 1
            info = Chem.AtomPDBResidueInfo(atomName=f'{element: >2}{counts[element]: <2}',
                                           residueName=self.resn,
                                           residueNumber=self.resi, chainId=self.chain)
            atom.SetPDBResidueInfo(info)

    def plonk(self, mol):
        """
        Temporarily here. Do not copy.
        There likely is a way to do this in OBabel
        This is using Fragmenstein ``MinimalPDBParser``.

        :param mol:
        :return:
        """
        pdbdata = MinimalPDBParser(self.pdb_block, remove_other_hetatms=True, ligname=self.resn)
        self.assign_pdb(mol)
        moldata = MinimalPDBParser(Chem.MolToPDBBlock(mol))
        pdbdata.append(moldata)
        return str(pdbdata)

    @singledispatchmethod
    def get_interaction_set(self) -> PLInteraction:
        """
        Overloaded method: block or mol return the iternaction set
        :return:
        """
        raise NotImplementedError

    @get_interaction_set.register
    def _(self, block: str) -> PLInteraction:
        holo = PDBComplex()
        holo.load_pdb(block, as_string=True)
        holo.analyze()
        return holo.interaction_sets[':'.join([self.resn, self.chain, str(self.resi)])]

    @get_interaction_set.register
    def _(self, mol: Chem.Mol) -> PLInteraction:
        if mol.GetNumAtoms() == 0:
            raise ValueError('Molecule has no atoms')
        holo = PDBComplex()
        holo.load_pdb(self.plonk(mol), as_string=True)
        holo.analyze()
        return holo.interaction_sets[':'.join([self.resn, self.chain, str(self.resi)])]

    def get_atomname(self, atom: Union[Atom, ob.OBAtom], atomnames: Optional[Sequence[str]]=None) -> str:
        """
        Given an atom, return its name.
        """
        if atomnames is not None and isinstance(atom, Atom):
            return atomnames[atom.idx - 1]  # Fortran indexing
        elif isinstance(atom, Atom):
            res: ob.OBResidue = atom.residue.OBResidue
            obatom = atom.OBAtom
        elif isinstance(atom, ob.OBAtom):
            obatom: ob.OBAtom = atom
            res: ob.OBResidue = obatom.GetResidue()
        else:
            raise TypeError
        return res.GetAtomID(obatom)  # this is likely to be ' C  ' as Babel has stripped this info

    def get_atom_by_atomname(self,
                             residue: Union[ob.OBResidue, Residue],
                             atomname: str,
                             atomnames: Optional[Sequence[str]]=None) -> ob.OBAtom:
        """
        Get an atom by its name in a residue.
        Note that the ligand will have its original atom names stripped by Babel.
        """
        if isinstance(residue, Residue):
            residue = residue.OBResidue
        obatom: ob.OBAtom
        for obatom in ob.OBResidueAtomIter(residue):
            if self.get_atomname(obatom, atomnames).strip() == atomname.strip():
                return obatom
        else:
            raise ValueError(f'No atom with name {atomname} in residue {residue.GetName()}')

    def get_interaction_counts(self, interaction_set: PLInteraction) -> Dict[Tuple[str, str, int], int]:
        """
        Count the number of interactions of each type for each residue
        """
        intxns: List = interaction_set.all_itypes
        intxn_dex = defaultdict(int)
        for intxn in intxns:
            key = (intxn.__class__.__name__, intxn.restype, intxn.resnr)
            intxn_dex[key] += 1  # noqa default dict works with tuples
        return dict(sorted(intxn_dex.items(), key=lambda kv: kv[0][2]))

    def summarize_interactions(self, atom_names: Sequence[str]) -> List[Dict[str, Any]]:
        interaction_set = self.get_interaction_set(self.pdb_block)
        details: List[Dict[str, Any]] = []
        for intxn in interaction_set.all_itypes:
            details.append(self.summarize_interaction(intxn))
        return details

    def summarize_interaction(self, intxn, atom_names: Sequence[str]) -> Dict[str, Any]:
        # https://github.com/openbabel/openbabel/blob/master/data/atomtyp.txt
        # https://blog.matteoferla.com/2023/07/a-note-no-plip-interactions.html
        relevant_atom_names = []
        type_name = intxn.__class__.__name__
        details = {'type': type_name, 'protein_resn': intxn.restype, 'protein_resi': intxn.resnr,
                   'protein_chain': intxn.reschain}
        if type_name == 'hbond':
            if intxn.protisdon:
                # assert hbond.a.residue.name != self.resn
                details['atom_names'] = [atom_names[intxn.a.idx - 1]]
                details['type'] = 'hbond_acceptor'
                details['babel_atom_types'] = [intxn.atype]
            else:
                details['atom_names'] = [atom_names[intxn.d.idx - 1]]
                details['type'] = 'hbond_donor'
                details['babel_atom_types'] = [intxn.dtype]
            details['distance'] = intxn.distance_ad  # intxn.distance_ad is to donor, _ah to hydrogen
        elif type_name == 'hydroph_interaction':
            details['atom_names'] = [atom_names[intxn.ligatom.idx - 1]]
            details['babel_atom_types'] = [intxn.ligatom.type]
            details['distance'] = intxn.distance
        elif type_name == 'pistack':
            details['atom_names'] = [atom_names[a.idx - 1] for a in intxn.ligandring.atoms]
            details['babel_atom_types'] = [a.type for a in intxn.ligandring.atoms]
            details['distance'] = intxn.distance
        elif type_name == 'waterbridge':
            if intxn.protisdon:
                # assert hbond.a.residue.name != self.resn
                details['atom_names'] = [atom_names[intxn.a.idx - 1]]
                details['type'] = 'water_acceptor'
                details['babel_atom_types'] = [intxn.atype]
                details['distance'] = intxn.distance_aw
            else:
                details['atom_names'] = [atom_names[intxn.d.idx - 1]]
                details['type'] = 'water_donor'
                details['babel_atom_types'] = [intxn.dtype]
                details['distance'] = intxn.distance_dw
        elif type_name == 'saltbridge':
            if intxn.protispos:
                details['atom_names'] = [atom_names[a.idx - 1] for a in intxn.negative.atoms]
                details['type'] = 'saltbridge_negative'
                details['babel_atom_types'] = [a.type for a in intxn.negative.atoms]
                details['distance'] = intxn.distance
            else:
                details['atom_names'] = [atom_names[a.idx - 1] for a in intxn.positive.atoms]
                details['type'] = 'saltbridge_positive'
                details['babel_atom_types'] = [a.type for a in intxn.positive.atoms]
                details['distance'] = intxn.distance
        elif type_name == 'pication':
            if intxn.protcharged:
                details['atom_names'] = [atom_names[a.idx - 1] for a in intxn.ring.atoms]
                details['type'] = 'pication_ring'
                details['babel_atom_types'] = [a.type for a in intxn.ring.atoms]
                details['distance'] = intxn.distance
            else:
                details['atom_names'] = [atom_names[a.idx - 1] for a in intxn.charge.atoms]
                details['type'] = 'pication_charge'
                details['babel_atom_types'] = [a.type for a in intxn.charge.atoms]
                details['distance'] = intxn.distance
        elif type_name == 'halogenbond':
            details['atom_names'] = [atom_names[intxn.don.idx - 1]]
            details['babel_atom_types'] = [intxn.don.type]
            details['distance'] = intxn.distance
        else:  # metal_complex basically.
            raise TypeError(type_name)
        return details