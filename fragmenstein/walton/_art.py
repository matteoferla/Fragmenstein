from __future__ import annotations
from io import StringIO
from typing import Dict, TYPE_CHECKING, Tuple, Union, List, Optional, Any
import sys

from ..display import MolNGLWidget, DISPLAYMODE, py3Dmol_monkey_patch, patched_3Dmol_view  # real or a mock
from IPython.display import display
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from ._base import WaltonBase


class WaltonArt(WaltonBase):

    def _mol2svg(self, mol: Chem.Mol) -> str:  # noqa it is static, but I dont care.
        d2d = Draw.rdMolDraw2D.MolDraw2DSVG(250, 200)
        mol2 = AllChem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol2)
        atoms = mol2.GetNumAtoms()
        for idx in range(atoms):
            mol2.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol2.GetAtomWithIdx(idx).GetIdx()))
        d2d.DrawMolecule(mol2)
        d2d.FinishDrawing()  # noqa it is filled.
        return d2d.GetDrawingText()  # noqa it is filled.

    def _get_mol_name(self, mol: Chem.Mol) -> str:  # noqa
        if mol.HasProp('_Name'):
            return mol.GetProp('_Name')
        else:
            return '「 no name 」'

    def _get_mol_color(self, mol: Chem.Mol) -> str:  # noqa
        if mol.HasProp('color'):
            return mol.GetProp('color')
        else:
            return 'gray'

    def _mol2smiles(self, mol: Chem.Mol) -> str:  # noqa
        mol2 = AllChem.RemoveHs(mol)
        return Chem.MolToSmiles(mol2)

    def _get_mol_details(self, mol: Chem.Mol) -> Dict[str, str]:
        return dict(name=self._get_mol_name(mol),
                    smiles=self._mol2smiles(mol),
                    svg=self._mol2svg(mol),
                    color=self._get_mol_color(mol),
                    )

    def _repr_html_(self):
        inner = '<p style="color: {color};"><b>{name}</b><br/>{smiles}</p>{svg}'
        template = f'<div style="float: left; padding: 10px;">{inner}</div>'
        return "\n".join(template.format(**self._get_mol_details(mol)) for mol in self.mols + [self.merged] if mol)

    # ## 3D --------------------------------------------

    def to_nglview(self) -> MolNGLWidget:
        """`
        generates a NGLWidget (``IPython.display.display`` will show it)
        with the compounds and the merged if present.
        The colours will be those in the ``Chem.Mol``'s property ``_color``.
        """
        self.monster.journal('This function is deprecated, use show instead')
        view = MolNGLWidget()
        for mol in self.mols + [self.merged]:
            if not mol:  # self.merged is empty
                continue
            view.add_mol(mol)
        return view

    def refresh_nglview(self, view: MolNGLWidget, *others:Chem.Mol) -> None:
        """
        In Walton altering hits manually does not have an effect on the NGLWidget,
        but altering hits has no bad juju unlike Monster where altering the `.mols` can have a detrimental effect,
        due to the lack of the sanitisation from the `monster.fix_mols(mols)` call.

        :param view:
        :return:
        """
        view.remove_all_components()
        for mol in self.mols + [self.merged] + list(others):
            if not mol:  # self.merged is empty
                continue
            view.add_mol(mol)


    def show3d(self) -> None:
        """
        Shows the structures both in 2d and 3d. Which is weird. Hence the deprecation.
        """
        self.monster.journal.warn('show3d is deprecated, use show instead')
        display(self)
        self.show()
        
    def show(self, to_display: bool = False, show_positioned_mol:bool=True, viewer_mode=DISPLAYMODE) \
                                                                -> Tuple[Union[MolNGLWidget,patched_3Dmol_view], str]:
        """
        Shows the structures in 3d and behaves like ``Monster.show``
        display(self)  will show in 2d molecules...
        """
        return self.monster.show(to_display=to_display, show_positioned_mol=show_positioned_mol, viewer_mode=viewer_mode)
