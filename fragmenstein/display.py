from typing import Sequence
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from IPython.display import display
from unittest.mock import Mock

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
    # monkey_patch is a function that needs to be called but it is attached as a cls method to py3Dmol.view
    # which for clarity is rebranded as ``patched_3Dmol_view``
    DISPLAYMODE = 'py3Dmol'
except Exception as error:
    patched_3Dmol_view = Mock(name='py3Dmol.view')


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
    hex_map = {atom.GetIdx(): atom.GetProp('_color') for atom in mol.GetAtoms()
               if atom.HasProp('_color') and atom.GetProp('_color')}
    rgb_map = {i: ColorConverter().to_rgb(n) for i, n in hex_map.items()}
    dopts.highlightAtomColors= list(rgb_map.values())
    dopts.highlightAtoms= list(rgb_map.keys())
    display(Draw.MolsToGridImage(flattos,
                                 legends=[mol.GetProp('_Name') if mol.HasProp('_Name') else '-' for mol in mols],
                                 subImgSize=subImgSize,
                                 useSVG=useSVG,
                                 molsPerRow=molsPerRow,
                                 drawOptions=dopts))

