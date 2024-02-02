# py3Dmol
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List

class py3Dmolplus(py3Dmol.view):
    # https://3dmol.csb.pitt.edu/doc/global.html#:~:text=as%20the%20colorscheme-,Properties,-Name
    color_names = ['default',
                     'greenCarbon',
                     'cyanCarbon',
                     'magentaCarbon',
                     'purpleCarbon',
                     'whiteCarbon',
                     'orangeCarbon',
                     'yellowCarbon',
                     'blueCarbon']
    def __init__(self, *args, **kwargs):
        """
        An extension of ``py3Dmol.view`` that adds control of the models by names.
        I.e. filling ``.model_addresses`` with the names of the models.
        The method ``.get_idx`` returns the index of the model with a given name.

        :param args:
        :param kwargs:
        """
        super().__init__(*args, **kwargs)
        self.model_addresses: List[str] = []

    def add_mol(self, mol, name=None, colorscheme='default', opacity=1., **kwargs):
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
        self.addModel(molblock, 'mol', name=name, **kwargs)
        # `.setStyle({'stick':{}})` will work, but I wanted the args to be explicit as I forget them.
        self.setStyle({'model': -1, }, {'stick': {'colorscheme': colorscheme, 'opacity': opacity}})

    def addModel(self, *args, name='unknown', **kwargs):
        super().addModel(*args, **kwargs)
        self.model_addresses.append(name)

    def get_idx(self, name):
        return self.model_addresses.index(name)
