__all__ = ['TestSet', 'MPro', 'BaseMolDataset']



from ._base_dataset_holder import BaseMolDataset, extend_doc
from rdkit import Chem
from typing import List
import os


# ----------------------------------------------------------------------------------------------------------------------
from . import test_mols

class TestSet(BaseMolDataset):
    """
    Test molecules for the tests.
    """
    dataset_package = test_mols

    @classmethod
    def get_template(cls):
        cls.raise_not_implemented()

    @classmethod
    def get_5SB7_mols(cls) -> List[Chem.Mol]:
        """
        In https://onlinelibrary.wiley.com/doi/10.1002/anie.202204052
        There is a lovely test case.
        Namely the followup bold4 is a merger of F36 and F04 hits.
        However... F04 (hello Hantzsch thiazole synthesis product) presents
        the problematic case that the thiazole can be flipped either way
        upon investigation of the crystal map
        and the followup bold4 is indeed a flipped merger.
        Therefore, without a custom_map the followup is incorrect.
        """
        return cls.get_sdf_mols('5SB7_mergers.sdf')

# ----------------------------------------------------------------------------------------------------------------------
from . import mpro_mols

class MPro(BaseMolDataset):
    """
    List of XChem hits of MPro from Fragalysis in Feb 2021.
    """

    dataset_package = mpro_mols

    @classmethod
    def file2mol_func(cls, filename: str) -> str:
        return os.path.splitext(filename)[0].replace('Mpro-', '').strip()

    @classmethod
    def mol2file_func(cls, hit_name: str) -> str:
        return f'Mpro-{hit_name}.mol'


# ----------------------------------------------------------------------------------------------------------------------
from . import mac1_mols

class Mac1(BaseMolDataset):
    """
    These are the fragments from https://pubmed.ncbi.nlm.nih.gov/33853786/
    binding he Macrodomain (Mac1 or NSP13) from SAR-COV-2
    """
    dataset_package = mac1_mols

extend_doc(TestSet)
extend_doc(MPro)
extend_doc(Mac1)



