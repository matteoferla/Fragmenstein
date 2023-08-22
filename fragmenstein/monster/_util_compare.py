from IPython.display import SVG, display, HTML
import nglview as nv
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from ..branding import divergent_colors
from functools import singledispatchmethod
from typing import Tuple, List, Dict, Optional
from io import StringIO
import re
from ..display import MolNGLWidget


def draw(mol, color_map, x=200, y=200, **kwargs):
    from matplotlib.colors import ColorConverter
    d = Draw.rdMolDraw2D.MolDraw2DSVG(x, y)
    d.drawOptions().addAtomIndices = True
    d.drawOptions().addStereoAnnotation = True
    d.drawOptions().prepareMolsBeforeDrawing = False
    d.drawOptions().dummiesAreAttachments = True
    flat = Chem.Mol(mol)
    AllChem.Compute2DCoords(flat)
    color_map2 = {i: ColorConverter().to_rgb(n) for i, n in color_map.items()}
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(d,
                                            flat,
                                            highlightAtoms=list(color_map.keys()),
                                            highlightAtomColors=color_map2,
                                            **kwargs
                                            )
    d.FinishDrawing()
    display(SVG(d.GetDrawingText()))



def get_idx(name, origin, default=None):
    rex = re.match(name + '\.(\d+)', origin)
    if rex:
        return int(rex.group(1))
    return default

class _MonsterUtilCompare:

    def get_color_origins(self) -> Dict[str, Dict[int, str]]:
        """
        Get for the hits and followup the color of the origin
        as seen in show_comparison
        :return:
        """
        color_maps = {}
        origins = self.origin_from_mol(self.positioned_mol)
        n_colors = sum(map(bool, origins))  # bool([]) is sensibly False in Python
        colors = divergent_colors[n_colors]
        mapped_idxs = [i for i, m in enumerate(origins) if m]
        followup_color_map = dict(zip(mapped_idxs, colors))
        positioned_name:str = self.positioned_mol.GetProp('_Name')
        color_maps[positioned_name] = followup_color_map
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
            color_maps[name] = color_map
        return color_maps

    def store_origin_colors_atomically(self):
        """
        Store the color of the origin in the mol as the private property _color
        :return:
        """
        # color_maps = {'benzene': {0: '#AED882', ...},..
        color_maps: Dict[str, Dict[int, str]] = self.get_color_origins()
        for hit in self.hits:  #:Chem.Mol
            name = hit.GetProp('_Name')
            for atom in hit.GetAtoms():  #: Chem.Atom
                atom.SetProp('_color', '#FFFFFF')  # white
            if name not in color_maps:
                continue
            for i in color_maps[name]:
                hit.GetAtomWithIdx(i % 100).SetProp('_color', color_maps[name][i])
        followup_map = color_maps[self.positioned_mol.GetProp('_Name')]
        for atom in self.positioned_mol.GetAtoms():  #: Chem.Atom
            atom.SetProp('_color', '#FFFFFF')  # white
        for i in followup_map:
            self.positioned_mol.GetAtomWithIdx(i).SetProp('_color', followup_map[i])


    def show_comparison(self, *args, **kwargs):
        """"
        Show the atom provenance of the follow molecule.

        Parameters:

        * nothing -> uses self.origin_from_mol
        * hit2followup -> uses the hit2followup dict
        """
        # ------- followup -----------------------
        #  origins is like [['MOL.7'], ['MOL.8'], ['MOL.9'], [], ...
        color_maps: Dict[str, Dict[int, str]] = self.get_color_origins()
        # -------------- hits ------------------
        for mol in self.hits:
            name = mol.GetProp('_Name')
            color_map = color_maps[name]
            print(f'hit {name}')  # legit print, not a debug scar
            draw(mol, color_map, 300, 300, **kwargs)
        print('Followup') # legit print, not a debug scar
        followup_color_map = color_maps[self.positioned_mol.GetProp('_Name')]
        draw(self.positioned_mol, followup_color_map, 300, 300)


    def convert_origins_to_custom_map(self, mol: Optional[Chem.Mol]=None, forbiddance=True) -> Dict[str, Dict[int, int]]:
        """
        The origins stored in the followup differ in format from the custom_map.
        The former is a list of lists of hit_name+atom_index,
        while the latter is a dictionary of hit_name to dictionary of
        hit atom indices to _intended_ followup index.
        This method converts the former to the latter.
        If forbiddance is True, non mapping atoms are marked with negatives.

        If mol is None, then self.positioned_mol is used.

        :return:
        """
        custom_map: Dict[str, Dict[int, int]] = {hit.GetProp('_Name'): {} for hit in self.hits}
        # `origins_from_mol` uses self.positioned_mol by default... what if it is not set?
        origins:List[List[str]] = self.origin_from_mol(mol)  # noqa It is in utility
        for fi, ori in enumerate(origins):
            for origin in ori:
                reg = re.search('(.*)\.(\d+)', str(origin))
                if reg is None:  # literally `none`... Although this is legacy.
                    continue
                name = reg.group(1)
                hi = int(reg.group(2)) % 100
                custom_map[name][hi] = fi
        new_custom_map = {}
        for hit in self.hits:
            name = hit.GetProp('_Name')
            ni = -1
            new_custom_map[name] = {}
            for hi in range(hit.GetNumAtoms()):
                if hi in custom_map[name]:
                    new_custom_map[name][hi] = custom_map[name][hi]
                elif not forbiddance:
                    pass
                elif hit.GetAtomWithIdx(hi).GetAtomicNum() == 1:
                    pass
                else:
                    new_custom_map[name][hi] = ni
                    ni -= 1
        return new_custom_map

    def _to_nglview_and_legend(self, show_positioned_mol:False) -> Tuple[nv.NGLWidget, str]:
        """
        This is called by both ``Monster.to_nglview`` and ``Victor.to_nglview``
        The color can be dictated by the optional private property ``_color``,
        which may have been assigned by ``walton.color_in()`` or the user.

        :return:
        """
        color_series = iter(divergent_colors[len(self.hits)])
        legend = ''
        view = MolNGLWidget()
        for mol in self.hits:
            colorValue = next(color_series) if not mol.HasProp('_color') else mol.GetProp('_color')
            view.add_mol(colorValue=colorValue, mol=mol)
            legend += f'<span style="color: {colorValue}">{mol.GetProp("_Name")}</span> '
        if show_positioned_mol and self.positioned_mol:
            view.add_mol(colorValue='white', mol=self.positioned_mol)
            legend += f'<span>positioned followup (in white)</span> '

        return view, legend

    def to_nglview(self, print_legend: bool = False) -> MolNGLWidget:
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

