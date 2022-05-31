from typing import Sequence
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from IPython.display import display
from matplotlib.colors import ColorConverter
import nglview as nv
from io import StringIO


def display_mols(mols: Sequence[Chem.Mol],
                 molsPerRow=5,
                 subImgSize=(150, 150),
                 useSVG=True):
    """
    Generic wrapper for calling ``display(Draw.MolsToGridImage)``
    """
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

class MolNGLWidget(nv.NGLWidget):
    """
    Adds a method ``add_mol`` that simply adds an rdkit molecule to the viewer.
    """

    def add_mol(self, mol: Chem.Mol, colorValue: str = '') -> nv.component.ComponentViewer:
        """
        Add a rdkit.Chem to an NGLWidget.
        This function is used by Walton and Monster

        :param mol: rdkit.Chem.Mol
        :param colorValue: if blank, the color in the property _color is used.
        :return:
        """
        if not mol:
            raise ValueError( 'Provided mol is None: '+\
                              'if user manual tinkering happened, please run monster.fix_hits')
        fh = StringIO(Chem.MolToPDBBlock(mol))  # I want atom names if present
        comp: nv.component.ComponentViewer = self.add_component(fh,  # noqa it's there.
                                                                name=mol.GetProp('_Name'),
                                                                ext='pdb'
                                                                )
        if not colorValue and mol.HasProp('_color'):
            colorValue = mol.GetProp('_color')
        if colorValue:
            comp.update_ball_and_stick(colorValue=colorValue)
        comp.update_ball_and_stick(multipleBond=True)
        # _color business stems from Walton.
        return comp

    def remove_all_components(self):
        self._js(f"""this.stage.removeAllComponents()""")
        # delattr() ???

    def add_neighbors(self, selection: str, comp_id: int = 0, radius: float = 5, style: str = 'hyperball',
                      color: str = 'gainsboro'):
        """
        Given a compounent id as interpreted by NGL JS, add the neighbours.

        :param selection:
        :param comp_id: if remove_all_components was called the component_0 etc are wrong.
        :param radius:
        :param style:
        :param color:
        :return:
        """
        self._js(f"""const comp = this.stage.compList[{comp_id}]
                     const target_sele = new NGL.Selection('{selection}');
                     const radius = {radius};
                     const neigh_atoms = comp.structure.getAtomSetWithinSelection( target_sele, radius );
                     const resi_atoms = comp.structure.getAtomSetWithinGroup( neigh_atoms );
                     comp.addRepresentation( "{style}", {{sele: resi_atoms.toSeleString(),
                     									  colorValue: "{color}",
                                                          multipleBond: true
                                                                       }});
                     comp.addRepresentation( "contact", {{sele: resi_atoms.toSeleString()}});
                """)
