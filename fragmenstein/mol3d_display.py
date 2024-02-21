"""
Due to the way ``py3Dmol`` works, I have to monkey patch it to add some functionality.
Namely, the ``py3Dmol.view.addModel`` is actually a __getattr__ call wrapping the JS,
which is not obeyed by super() and thus cannot be overridden by subclassing.

Therefore, I have to monkey patch it.

... code-block:: python

    viewer = py3Dmol.view()
    viewer.monkey_patch()
    viewer.add_template(vicky.apo_pdbblock)
    viewer.add_mol(vicky.hits[0], name='hit', carbon_color='cyan')
    viewer.add_mol(vicky.minimized_mol, name=vicky.long_name, colorscheme='coralCarbon')
    viewer.show()

carbon_color can be None (default), a CSS color name, or a hex color.
"""

# py3Dmol
import py3Dmol
from types import MethodType
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List, Optional

# https://3dmol.csb.pitt.edu/doc/global.html#:~:text=as%20the%20colorscheme-,Properties,-Name
py3Dmol.view.color_names = ['default',
                            'greenCarbon',
                            'cyanCarbon',
                            'magentaCarbon',
                            'purpleCarbon',
                            'whiteCarbon',
                            'orangeCarbon',
                            'yellowCarbon',
                            'blueCarbon']


def _add_mol(self, mol, name=None, carbon_color: Optional[str] = None, opacity=1., **kwargs):
    """
    Whereas there is ``IPythonConsole.addMolToView(caffeine, view)`` in rdkit using py3Dmol,
    I want more.

    :param mol:
    :param name:
    :param colorscheme:
    :param opacity:
    :param kwargs:
    :return:
    """
    if name is not None:
        pass
    elif mol.HasProp('_Name'):
        name: str = mol.GetProp('_Name')
    else:
        raise ValueError('No name provided and mol has no _Name property (title)')
    if name in self.model_addresses:
        raise ValueError(f'{name} is already present in the viewer')
    molblock: str = Chem.MolToMolBlock(mol)
    self.add_model(molblock, 'mol', name=name, **kwargs)
    # `.setStyle({'stick':{}})` will work, but I wanted the args to be explicit as I forget them.
    if carbon_color == 'default':
        self.setStyle({'model': -1, }, {'stick': {'colorscheme': 'default', 'opacity': opacity}})
    elif not isinstance(carbon_color, str):
        raise ValueError(f'No idea what is carbon_color={carbon_color}')
    elif 'Carbon' in carbon_color:
        self.setStyle({'model': -1, }, {'stick': {'colorscheme': carbon_color, 'opacity': opacity}})
    elif '#' not in carbon_color:
        self.setStyle({'model': -1, }, {'stick': {'colorscheme': f'{carbon_color}Carbon', 'opacity': opacity}})
    elif carbon_color == 'prop' or carbon_color is None:
        color = mol.GetProp('color')
        self.setStyle({'model': -1, }, {'stick': {'colorscheme': f'{color}Carbon', 'opacity': opacity}})
    else:
        # this will not work on an update which runs off updatejs
        self.startjs += f'''\n
    let customColorizeFor{name} = function(atom){{
          // attribute elem is from https://3dmol.csb.pitt.edu/doc/AtomSpec.html#elem
          if (atom.elem === 'C'){{
              return "{carbon_color}"
          }}else{{
              return $3Dmol.getColorFromStyle(atom, {{colorscheme: "whiteCarbon"}});
          }}
      }}\n'''
        self.setStyle({'model': -1}, {'stick': {'colorfunc': 'customColorize', 'opacity': opacity}})
        # make it a function not a string "customColorize"
        self.startjs = self.startjs.replace('"customColorize"', f'customColorizeFor{name}')


def _add_template(self, pdbblock: str, name='template', colorscheme='whiteCarbon'):
    self.add_model(pdbblock, 'PDB', name=name)
    self.setStyle({'model': -1}, {'cartoon': {'colorscheme': colorscheme, 'opacity': 1}})
    self.zoomTo()


def _get_idx(self, name):
    return self.model_addresses.index(name)


def _add_model(self, *args, name='unknown', **kwargs):
    self.addModel(*args, **kwargs)
    self.model_addresses.append(name)


def monkey_patch(viewer: py3Dmol.view):
    viewer.model_addresses: List[str] = []
    viewer.add_model = _add_model.__get__(viewer, py3Dmol.view)
    viewer.add_mol = _add_mol.__get__(viewer, py3Dmol.view)
    viewer.get_idx = _get_idx.__get__(viewer, py3Dmol.view)
    viewer.add_template = _add_template.__get__(viewer, py3Dmol.view)

py3Dmol.view.monkey_patch = monkey_patch

patched_3Dmol_view = py3Dmol.view