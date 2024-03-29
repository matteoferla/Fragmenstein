from typing import (Union, Dict, List)

from rdkit import Chem
from rdkit.Chem import AllChem

from ..monster import Monster
from ..branding import divergent_colors
from ..display import color_in

# courtesy of .legacy
from functools import singledispatchmethod

class WaltonBase:
    color_scales = divergent_colors  # NO LONGER USABLE TODO: remove

    def __init__(self,
                 mols: List[Chem.Mol],
                 superposed: bool = False):
        """
        To initialised from SMILES use the classmethod ``.from_smiles``.
        These are assumed to have a conformer.
        The mols will be assigned a property ``_color`` based on
        the class attribute color_scales. By default it uses ``fragmenstein.branding.divergent_colors``

        :param mols: list of mols
        :param superposed: are they superposed? sets the namesake argument that does nothing ATM
        """
        # ## Mol
        self.mols: List[Chem.Mol] = mols  # It does not alter then, but `monster.fix_hits` will be called on __call__...
        # assign index (just in case):
        for idx, mol in enumerate(self.mols):
            mol.SetIntProp('_mol_index', idx)
        # assign color:
        color_in(self.mols, skip_feija=True)
        # ## superposed
        self.superposed: bool = superposed
        # ## Computed
        self.merged: Union[None, Chem.Mol] = None  # not a dynamic property: allows vandalism
        # the hits will be passed again on call:
        self.monster = Monster(list(map(AllChem.RemoveHs, self.mols)))

    def color_in(self, skip_feija=False):
        """
        assigns a color property to a mol based on color_scales of correct length

        The first colour is the Fragmenstein colour (feijoa). Setting `color_in(False)` will skip it,
        allowing it to be used later on.

        Gets called by ``__init__`` and ``duplicate``.
        """
        self.monster.journal.warning('color_in as a Walton method is deprecated. Use the global color_in')
        color_in(self.mols, skip_feija=skip_feija)

    @classmethod
    def from_smiles(cls,
                    superposed: bool = False,
                    add_Hs: bool = False,
                    **name2smiles: str):
        """
        Load from SMILES.
        provided as named arguments: ``from_smiles(bezene='c1ccccc1',..)``
        """
        mols: List[Chem.Mol] = []
        for name, smiles in name2smiles.items():
            mol = Chem.MolFromSmiles(smiles)
            if add_Hs:
                mol = AllChem.AddHs(mol, addCoords=True)
            AllChem.EmbedMolecule(mol)
            mol.SetProp('_Name', name)
            mols.append(mol)
        return cls(mols=mols, superposed=superposed)

    def __call__(self, color='#a9a9a9', minimize: bool = False, **combine_kwargs) -> Chem.Mol:  # darkgrey
        """
        Calls Monster to do the merger.
        Filling the attribute ``merged`` w/ a Chem.Mol.
        Also returns it.
        """
        # neogreen '#39ff14'
        # joining_cutoff= 5
        self.monster.hits = list(map(AllChem.RemoveHs, self.mols))
        self.monster.combine(**combine_kwargs)
        if minimize:
            self.monster.positioned_mol = self.monster.mmff_minimize(allow_lax=True).mol
        self.monster.store_origin_colors_atomically()
        self.merged = AllChem.RemoveHs(self.monster.positioned_mol)
        self.merged.SetProp('color', color)
        return self.merged

    def duplicate(self, mol_idx: int):
        """
        Duplicate the molecule at a given index.
        And fix colours.
        """
        self.mols.append(Chem.Mol(self.get_mol(mol_idx)))
        color_in(self.mols, skip_feija=True)

    @singledispatchmethod
    def get_mol(self, mol_idx: int) -> Chem.Mol:
        """
        Type dispatched method:

        * Gets the molecule in ``.mols`` with index ``mol_idx``
        * returns the molecule provided as ``mol``

        The latter route is not used within the module
        but does mean one could pass a mol instead of a mol_idx...
        """
        assert isinstance(mol_idx, int)
        assert len(self.mols) > mol_idx, f'The instance of Walton has only {len(self.mols)}, so cannot get {mol_idx}.'
        return self.mols[mol_idx]

    @get_mol.register
    def _(self, mol: Chem.Mol) -> Chem.Mol:
        return mol
