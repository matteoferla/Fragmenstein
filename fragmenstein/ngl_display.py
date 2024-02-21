"""
NGL may not be installed: cf https://github.com/matteoferla/Fragmenstein/issues/37
as it's an optional: NGLview breaks Colaboratory in Jan 2024.

Both Victor, Monster and Walton use this optionally via the method ``.to_nglview``
which is inherited from:
_MonsterUtilCompare.to_nglview
_VictorShow.to_nglview
WaltonArt.to_nglview
"""

from rdkit import Chem
import nglview as nv
from nglview.component import ComponentViewer
from io import StringIO
import warnings

class MolNGLWidget(nv.NGLWidget):
    """
    Adds:

    a method ``add_mol`` that simply adds an rdkit molecule to the viewer.

    .. code-block:: python

        from rdkit import Chem
        mol = Chem.MolFromSmiles('CC')
        view = MolNGLWidget()
        view.add_mol(mol, colorValue='pink')

    A method ``add_neighbors`` that adds a representation of the neighbours of a given selection.

    .. code-block:: python

        view.add_neighbors('[ATP]', comp_id=0, radius=5, style='hyperball', color='gainsboro')

    """

    def add_mol(self, mol: Chem.Mol, carbon_color: str = '', colorValue:str='') -> ComponentViewer:
        """
        Add a rdkit.Chem to an NGLWidget.
        This function is used by Walton and Monster

        :param mol: rdkit.Chem.Mol
        :param colorValue: if blank, the color in the property _color is used.
        :return:
        """
        if carbon_color:
            # this is to standardise the color names with py3Dmol
            colorValue = carbon_color
        if not mol:
            raise ValueError('Provided mol is None: ' + \
                             'if user manual tinkering happened, please run monster.fix_hits')
        fh = StringIO(Chem.MolToPDBBlock(mol))  # I want atom names if present
        comp: ComponentViewer = self.add_component(fh,  # noqa it's there.
                                                   name=mol.GetProp('_Name'),
                                                   ext='pdb'
                                                   )
        comp.remove_ball_and_stick()
        if not colorValue and mol.HasProp('color'):
            colorValue = mol.GetProp('color')
        if colorValue:
            comp.add_representation('ball+stick', colorValue=colorValue, multipleBond=True)
        else:
            comp.add_representation('ball+stick', multipleBond=True)
        # _color business stems from Walton.
        return comp

    def remove_all_components(self):
        """
        This removes the components and all traces of their existance in the widget

        :return:
        """
        self._js(f"""this.stage.removeAllComponents()""")
        i = 0
        while hasattr(self, f'component_{i}'):
            delattr(self, f'component_{i}')
            i += 1
        self._ngl_component_ids = []

    def add_neighbors(self, selection: str, comp_id: int = 0, radius: float = 5, style: str = 'hyperball',
                      color: str = 'gainsboro'):
        """
        Given a compounent id as an integer interpreted by NGL JS (not 'component_0'), show the neighbours.

        :param selection: NGL style selection string, e.g. ``[ATP] or 1-30:A``
        :param comp_id: integer not string
        :param radius: float of atom to atom distance max
        :param style: see NGL manual for the representation styles. Default is ``hyperball``
        :param color: hashtag prefixed hex color string or CSS color name
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

    def add_selection_signal(self, molname_id: str, atom_id: str):
        """
        Add a signal to the viewer than fills the elements #molname_id and #atom_id
        when user clicks on an atom —not bond.


        .. code-block:: python

            from IPython.display import display, HTML
            # in real code the html entities &gt;/&lt; would be actual greater-than and lesser-than signs
            display(HTML(f'&gt;div id="#{molname_id}"&lt;&gt;/div>&gt;div id="#{atom_id}"&lt;&gt;/div&lt;'))
            view.add_selection_signal(molname_id, atom_id)

        :param molname_id:
        :param atom_id:
        :return:
        """
        self._js(f'''this.stage.signals.clicked.add(pickingProxy => {{
                if (pickingProxy && (pickingProxy.atom || pickingProxy.bond )){{
                    const atom = pickingProxy.atom || pickingProxy.closestBondAtom;
                    const component = pickingProxy.component;
                    console.log(atom);
                    console.log(component);
                    document.getElementById('{molname_id}').innerText = `${{component.name}}`;
                    document.getElementById('{atom_id}').innerText = `${{atom.index}}`;
                    }}
            }});
        ''')