import os
from io import StringIO, BytesIO
from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import List, Optional
from types import ModuleType
import importlib.resources as pkg_resources
import random


class BaseMolDataset:
    """
    Helper functions for loading the data.
    I had considered making it a "module factory" but gets quickly unintelligible.
    """

    dataset_package: ModuleType = '__main__'

    def __init__(self):
        raise NotImplementedError('Virtual method')

    def raise_not_implemented(self):
        raise NotImplementedError('Method not implemented for this dataset')

    @classmethod
    def file2mol_func(cls, filename: str) -> str:
        return os.path.splitext(filename)[0].strip()

    @classmethod
    def mol2file_func(cls, hit_name: str) -> str:
        return f'{hit_name.strip()}.mol'

    @classmethod
    def get_molblock(cls, hit_name: Optional[str] = None) -> str:
        """
        returns the mol block for the given hit. Chosen randomly if unspecified.
        See ``get_mpro_hit_list()`` for options.
        """
        if hit_name is None:
            hit_name = random.choice(cls.get_hit_list())
        return cls.get_text(cls.mol2file_func(hit_name))
    
    @classmethod
    def get_hit_list(cls) -> List[str]:
        return [cls.file2mol_func(fn) for fn in pkg_resources.contents(cls.dataset_package) if '.mol' in fn]
    
    @classmethod
    def get_text(cls, filename: str) -> str:
        return pkg_resources.read_text(cls.dataset_package, filename)

    @classmethod
    def get_bytes(cls, filename: str) -> bytes:
        return pkg_resources.read_binary(cls.dataset_package, filename)
    
    @classmethod
    def exists(cls, filename: str) -> bool:
        return pkg_resources.is_resource(cls.dataset_package, filename)

    @classmethod
    def get_mol(cls, hit_name: Optional[str] = None) -> Chem.Mol:
        """
        returns the Chem.Mol instance for the given hit. Chosen randomly if unspecified.
        See ``get_mpro_hit_list()`` for options.
        """
        mol = Chem.MolFromMolBlock(cls.get_molblock(hit_name))
        mol.SetProp('_Name', hit_name)
        return mol

    @classmethod
    def get_sdf_mols(cls, filename: str) -> List[Chem.Mol]:
        """
        Get the Chem.Mols in the specified sdf file.

        :param filename:
        :return:
        """
        bytes_io = BytesIO(cls.get_bytes(filename))
        with Chem.ForwardSDMolSupplier(bytes_io) as reader:
            return list(reader)

    @classmethod
    def get_filtered_mol(cls, mw_cutoff: float = float('inf')) -> Chem.Mol:
        """
        Return a random Chem.Mol from the MPro dataset,
        with one of the following filtering criteria:

        * ``mw_cutoff``: molWt cutoff
        """
        # For now only MolWt
        choices = cls.get_hit_list()
        random.shuffle(choices)
        for hit_name in choices:
            mol = cls.get_mol(hit_name)
            if mw_cutoff >= Descriptors.MolWt(mol):
                return mol
        else:
            raise ValueError('No hit matches the specified values')

    @classmethod
    def get_n_filtered_mols(cls, amount: int, **cutoffs) -> List[Chem.Mol]:
        """Get ``amount`` of the mols (Chem.Mol) randomly
         that match a cutoff criterion.
        As listed in ``get_filtered_mol``"""
        mols = []
        mol_names = set()
        assert amount > 0, 'A zero amount does nothing.'
        for i in range(len(cls.get_hit_list())):
            mol = cls.get_filtered_mol(**cutoffs)
            mol_name = mol.GetProp('_Name')
            if mol_name in mol_names:
                continue
            mol_names.add(mol_name)
            mols.append(mol)
            if len(mols) == int(amount):
                return mols
        else:
            raise ValueError(f'There are only {len(mols)} molecules ({mol_names}) that match the criteria {cutoffs}')

    @classmethod
    def get_template(cls) -> str:
        return cls.get_text('template.pdb')

# not for the classes to use
def extend_doc(cls):
    """
    This adds usage notes dynamically to the modules.
    The spacing will be off.
    """
    extra = f"""
    The module {cls.dataset_package.__name__} contains the data, but
    not functions to retrieve them without using ``pkg_resources``.
    These are added by the class {cls.__name__} as class methods.
    Not all of the methods are implemented for all datasets.
    """
    cls.dataset_package.__doc__ += extra
