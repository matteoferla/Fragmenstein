from IPython.display import SVG, display
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from matplotlib.colors import ColorConverter
from ..branding import divergent_colors
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


def get_idx(name, origin):
    rex = re.match(name + '\.(\d+)', origin)
    if rex:
        return int(rex.group(1))

class _MonsterUtilCompare:

    def show_comparison(self):
        """"
        Show the atom provenance of the follow molecule
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
            print(f'hit {name}')
            draw(mol, color_map, 300, 300)
        print('Followup')
        draw(self.positioned_mol, followup_color_map, 300, 300)
