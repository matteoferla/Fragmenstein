from IPython.display import SVG, display, HTML
import nglview as nv
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from matplotlib.colors import ColorConverter
from ..branding import divergent_colors
from functools import singledispatchmethod
from typing import Tuple, List, Dict
from io import StringIO
import re


def draw(mol, color_map, x=200, y=200):
    d = Draw.rdMolDraw2D.MolDraw2DSVG(x, y)
    flat = Chem.Mol(mol)
    AllChem.Compute2DCoords(flat)
    color_map2 = {i: ColorConverter().to_rgb(n) for i, n in color_map.items()}
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(d,
                                            flat,
                                            highlightAtoms=list(color_map.keys()),
                                            highlightAtomColors=color_map2,
                                            )
    d.FinishDrawing()
    display(SVG(d.GetDrawingText()))


def get_idx(name, origin, default=None):
    rex = re.match(name + '\.(\d+)', origin)
    if rex:
        return int(rex.group(1))
    return default

class _MonsterUtilCompare:

    def show_comparison(self):
        """"
        Show the atom provenance of the follow molecule.

        Parameters:

        * nothing -> uses self.origin_from_mol
        * hit2followup -> uses the hit2followup dict
        """
        # ------- followup -----------------------
        #  origins is like [['MOL.7'], ['MOL.8'], ['MOL.9'], [], ...
        origins = self.origin_from_mol(self.positioned_mol)
        n_colors = sum(map(bool, origins))  # bool([]) is sensibly False in Python
        colors = divergent_colors[n_colors]
        mapped_idxs = [i for i, m in enumerate(origins) if m]
        # d = Draw.rdMolDraw2D.MolDraw2DCairo(500, 500)
        followup_color_map = dict(zip(mapped_idxs, colors))
        # later: draw(victor.minimized_mol, color_map)

        # -------------- hits ------------------
        for mol in self.hits:
            name = mol.GetProp('_Name')
            # convert ['MOL.7'] to a list of indices of MOL
            matched = [[get_idx(name, o) for o in ori] for ori in origins]
            # make a color_map
            color_map = {}
            for followup_idx, target_idxs in enumerate(matched):
                for i in target_idxs:  # there ought to only zero or one...
                    if i is None:
                        continue
                    color_map[i] = followup_color_map[followup_idx]
            print(f'hit {name}')  # legit print, not a debug scar
            draw(mol, color_map, 300, 300)
        print('Followup') # legit print, not a debug scar
        draw(self.positioned_mol, followup_color_map, 300, 300)


    def convert_origins_to_custom_map(self) -> Dict[str, Dict[int, int]]:
        """
        The origins stored in the followup differ in format from the custom_map.
        The former is a list of lists of hit_name+atom_index,
        while the latter is a dictionary of hit_name to dictionary of
        hit atom indices to _intended_ followup index.
        This method converts the former to the latter.

        :return:
        """
        custom_map: Dict[str, Dict[int, int]] = {}
        # `origins_from_mol` uses self.positioned_mol by default... what if it is not set?
        origins:List[List[str]] = self.origin_from_mol()  # noqa It is in utility
        for hit in self.hits:  #:Chem.Mol
            name:str = hit.GetProp('_Name')
            # "default=-1-fi" is to assign a unique negative number to prevent anything mapping to the target index.
            custom_map[name] = {get_idx(name, o, default=-1-fi): fi for fi, ori in enumerate(origins) for o in ori}
        return custom_map

    def _to_nglview_and_legend(self, show_positioned_mol:False) -> Tuple[nv.NGLWidget, str]:
        """
        This is called by both ``Monster.to_nglview`` and ``Victor.to_nglview``
        :return:
        """
        color_series = iter(divergent_colors[len(self.hits)])
        legend = ''
        view = nv.NGLWidget()
        for mol in self.hits:
            colorValue = next(color_series) if not mol.HasProp('_color') else mol.GetProp('_color')
            self.add_mol_to_nglview(colorValue=colorValue, mol=mol, view=view)

            legend += f'<span style="color: {colorValue}">{mol.GetProp("_Name")}</span> '
        if show_positioned_mol and self.positioned_mol:
            self.add_mol_to_nglview(colorValue='white', mol=self.positioned_mol, view=view)
            legend += f'<span>positioned followup (in white)</span> '

        return view, legend

    def add_mol_to_nglview(self, view: nv.NGLWidget, mol: Chem.Mol, colorValue: str) -> nv.component.ComponentViewer:
        if not mol:
            raise ValueError(
                'One of the hits is None: if user manual tinkering happened, please run monster.fix_hits')
        fh = StringIO(Chem.MolToMolBlock(mol))  # I want atom names
        comp: nv.component.ComponentViewer = view.add_component(fh,  # noqa it's there.
                                                                name=mol.GetProp('_Name'),
                                                                ext='mol')
        # _color business stems from Walton.
        comp.update_ball_and_stick(colorValue=colorValue, multipleBond=True)
        return comp

    def to_nglview(self, print_legend: bool = False) -> nv.NGLWidget:
        """
        This is not the same method as in Victor.
        generates a NGLWidget (``IPython.display.display`` will show it)
        with the compounds and the merged if present.
        To override the colours:
        The colours will be those in the ``Chem.Mol``'s property ``_color`` if present.

        Returns -> nv.NGLWidget
        """
        view, legend = self._to_nglview_and_legend(show_positioned_mol=True)
        if print_legend:
            display(HTML(legend))
        # async madness: disabled for now.
        # view.center(f'[{self.ligand_resn}]')
        return view

