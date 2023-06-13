## This is an optional import

"""
Trimmed down version of PLIPper from PLIP-PyRosetta-hotspots-test
https://github.com/matteoferla/PLIP-PyRosetta-hotspots-test/blob/main/plipspots_docking/plipspots/serial.py
"""

from collections import defaultdict
from functools import cached_property
from typing import Tuple, Dict, List, Union

from openbabel.pybel import Atom, Residue
from openbabel.pybel import ob
from plip.structure.preparation import PDBComplex, PLInteraction


class PLIPper:
    def __init__(self, pdb_block: str, resn: str, resi: int, chain: str):
        assert 'ATOM' in pdb_block, f'No ATOM entry in block provided: {pdb_block}'
        self.pdb_block = pdb_block # victor.minimized_pdbblock
        self.resn = resn  # victor.ligand_resn
        self.resi = resi  # victor.ligand_resi[:-1]
        self.chain = chain  # victor.ligand_resi[-1]

    @cached_property
    def interaction_set(self) -> PLInteraction:
        holo = PDBComplex()
        holo.load_pdb(self.pdb_block, as_string=True)
        holo.analyze()
        return holo.interaction_sets[(self.resn, self.chain, str(self.resi))]

    @cached_property
    def interaction_counts(self) -> Dict[Tuple[str, str, int], int]:
        """
        Count the number of interactions of each type for each residue
        """
        intxns: List = self.interaction_set.all_itypes
        intxn_dex: Dict[Tuple[str, str, int], int] = defaultdict(int)
        for intxn in intxns:
            key = (str(intxn.__class__.__name__), str(intxn.restype), int(intxn.resnr))
            intxn_dex[key] += 1
        return dict(sorted(intxn_dex.items(), key=lambda kv: kv[0][2]))

    # these are not needed... but can be useful
    def get_atomname(self, atom: Union[Atom, ob.OBAtom]) -> str:
        """
        Given an atom, return its name.
        """
        if isinstance(atom, Atom):
            res: ob.OBResidue = atom.residue.OBResidue
            obatom = atom.OBAtom
        elif isinstance(atom, ob.OBAtom):
            obatom: ob.OBAtom = atom
            res: ob.OBResidue = obatom.GetResidue()
        else:
            raise TypeError
        return res.GetAtomID(obatom)

    def get_atom_by_atomname(self, residue: Union[ob.OBResidue, Residue], atomname: str) -> ob.OBAtom:
        """
        Get an atom by its name in a residue.
        """
        if isinstance(residue, Residue):
            residue = residue.OBResidue
        obatom: ob.OBAtom
        for obatom in ob.OBResidueAtomIter(residue):
            if residue.GetAtomID(obatom).strip() == atomname:
                return obatom
        else:
            raise ValueError(f'No atom with name {atomname} in residue {residue.GetName()}')
