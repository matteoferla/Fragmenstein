from typing import Sequence
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from IPython.display import display
from unittest.mock import Mock
from typing import List

DISPLAYMODE = 'rdkit'
try:
    # Optional. See docstring in .ngl_display for more.
    from .ngl_display import MolNGLWidget, nv, ComponentViewer
    NGLWidget = nv.NGLWidget
    DISPLAYMODE = 'ngl'
except Exception as error:
    MolNGLWidget = Mock(name='MolNGLWidget')
    ComponentViewer = Mock(name='nglview.component.ComponentViewer')
    NGLWidget = Mock(name='nglview.NGLWidget')

try:
    from .mol3d_display import patched_3Dmol_view
    from .mol3d_display import monkey_patch as py3Dmol_monkey_patch
    # monkey_patch is a function that needs to be called but it is attached as a cls method to py3Dmol.view
    # which for clarity is rebranded as ``patched_3Dmol_view``
    DISPLAYMODE = 'py3Dmol'
except Exception as error:
    patched_3Dmol_view = Mock(name='py3Dmol.view')
    py3Dmol_monkey_patch = Mock(name='py3Dmol.view.monkey_patch')
from .branding import divergent_colors


def color_in(mols: List[Chem.Mol], color_scale=None, skip_feija=False):
    """
    assigns a color property to a mol based on color_scales of correct length

    In the `divergent_colors` first colour is the Fragmenstein colour (feijoa). Setting `color_in(False)` will skip it,
    allowing it to be used later on.
    """
    n_mols = len(mols)
    if n_mols == 0:
        return
    if color_scale is None and skip_feija:
        color_scale = divergent_colors[n_mols + int(skip_feija)][1:]
    elif color_scale is None:
        color_scale = divergent_colors[n_mols]
    elif len(color_scale) < n_mols:
        raise ValueError(f'color_scale is too short for {n_mols} mols.')
    else:
        pass
    for mol, color in zip(mols, color_scale):
        mol.SetProp('color', color)


def display_mols(mols: Sequence[Chem.Mol],
                 molsPerRow=5,
                 subImgSize=(150, 150),
                 useSVG=True):
    """
    Generic wrapper for calling ``display(Draw.MolsToGridImage)``
    """
    from matplotlib.colors import ColorConverter

    if not mols:
        return  # no mols, no display
    flattos = [AllChem.RemoveHs(mol) for mol in mols if isinstance(mol, Chem.Mol)]
    for mol in flattos:
        AllChem.Compute2DCoords(mol)
    dopts = Draw.rdMolDraw2D.MolDrawOptions()  # noqa
    dopts.addAtomIndices = True
    hex_map = {atom.GetIdx(): atom.GetProp('color') for atom in mol.GetAtoms()
               if atom.HasProp('color') and atom.GetProp('color')}
    rgb_map = {i: ColorConverter().to_rgb(n) for i, n in hex_map.items()}
    dopts.highlightAtomColors= list(rgb_map.values())
    dopts.highlightAtoms= list(rgb_map.keys())
    display(Draw.MolsToGridImage(flattos,
                                 legends=[mol.GetProp('_Name') if mol.HasProp('_Name') else '-' for mol in mols],
                                 subImgSize=subImgSize,
                                 useSVG=useSVG,
                                 molsPerRow=molsPerRow,
                                 drawOptions=dopts))

