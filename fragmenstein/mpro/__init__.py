########################################################################################################################

__doc__ = \
    """
This is a variant of Victor for MPro that uses data from PostEra
    """

########################################################################################################################


from ..victor import Victor
from . import data

import os, pyrosetta
from rdkit import Chem
from typing import List, Optional

import pandas as pd
import io
import requests
import random
import importlib.resources as pkg_resources

def get_mpro_template():
    return pkg_resources.read_text(data, 'template.pdb')

def _clean_hitname(hit_name: Optional[str]=None) -> str:
    if hit_name is None:
        hit_name = random.choice(get_mpro_hit_list())
        'Mpro-6y2g.mol'
    hit_name = hit_name.replace('.mol', '').replace('Mpro-', '').strip()
    if not pkg_resources.is_resource(data.hit_mols, f'Mpro-{hit_name}.mol'):
        raise FileNotFoundError(f'There is no hit {hit_name} in the cache. Choices are {get_mpro_hit_list()}')
    return hit_name

def get_mpro_molblock(hit_name: Optional[str]=None) -> str:
    """
    returns the mol block for the given hit. Chosen randomly if unspecified.
    See ``get_mpro_hit_list()`` for options.
    """
    hit_name = _clean_hitname(hit_name)
    return pkg_resources.read_text(data.hit_mols, f'Mpro-{hit_name}.mol')

def get_mpro_mol(hit_name: Optional[str]=None) -> Chem.Mol:
    """
    returns the Chem.Mol instance for the given hit. Chosen randomly if unspecified.
    See ``get_mpro_hit_list()`` for options.
    """
    hit_name = _clean_hitname(hit_name)
    mol = Chem.MolFromMolBlock(get_mpro_molblock(hit_name))
    mol.SetProp('_Name', hit_name)
    return mol

def get_mpro_hit_list():
    """
    List of XChem hits of MPro from Fragalysis in Feb 2021.
    """
    return [fn.replace('.mol', '').replace('Mpro-', '') for fn in pkg_resources.contents(data.hit_mols) if '.mol' in fn]


def pose_fx(pose: pyrosetta.Pose):
    """
    Histidine in delta.
    """
    pdb2pose = pose.pdb_info().pdb2pose
    r = pdb2pose(res=41, chain='A')
    MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
    MutateResidue(target=r, new_res='HIS').apply(pose)


def poised_pose_fx(pose: pyrosetta.Pose):
    """
    Histidine in delta and cysteine in thiolate.
    """
    pdb2pose = pose.pdb_info().pdb2pose
    r = pdb2pose(res=41, chain='A')
    MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
    MutateResidue(target=r, new_res='HIS_D').apply(pose)
    r = pdb2pose(res=145, chain='A')
    MutateResidue(target=r, new_res='CYZ').apply(pose)

class MProVictor(Victor):
    constraint_function_type = 'FLAT_HARMONIC'

    @classmethod
    def from_hit_codes(cls, hit_codes: List[str], **options):
        hits = [get_mpro_mol(xnumber) for xnumber in hit_codes]
        return cls(hits=hits, **options)

    def __init__(self, category:Optional[str]=None, **options):
        # this category flag is solely for Mpro?
        # it stems from the Moonshot file.
        self.category = category
        if category == 'noncolavent':
            fx = poised_pose_fx
        else:
            fx = pose_fx
        # --------------------
        defaults = dict(
            pdb_block=get_mpro_template(),
            ligand_resn='LIG',
            ligand_resi='1B',
            covalent_resn='CYS',
            covalent_resi='145A',
            extra_protein_constraint='AtomPair  SG  145A  NE2  41A HARMONIC 3.5 0.2\n',
            pose_fx=fx  # from the above.
        )
        super().__init__(**{**defaults, **options})

    def _determine_extras(self, smiles):
        defaults = dict()
        if self.category not in (None, 'noncovalent') and '_' in self.category:
            cname, rxd = self.category.split('_')
            if rxd == 'noncovalent':
                wd = [wd for wd in self.warhead_definitions if wd['name'] == cname][0]
                mol = Chem.MolFromSmiles(smiles)
                nc = Chem.MolFromSmiles(wd['noncovalent'])
                atomnames = dict(zip(mol.GetSubstructMatch(nc), wd['noncovalent_atomnames']))
                extra_constraint = 'AtomPair  SG  145A  CX   1B HARMONIC 3.2 0.5\n'
                extra_constraint += wd['constraint']
                defaults = dict(atomnames=atomnames,
                                extra_ligand_constraint=extra_constraint)
        return defaults

    #self.combine(**options) unchanged.

    def place(self, **options):
        defaults = self._determine_extras(options['smiles'])
        return super().place(**{**defaults, **options})

    # ======= postera csv file ops =====================================================================================

    @classmethod
    def fetch_postera(cls):
        """
        Reads the submission file off Github.
        For a local version, just ``postera = pd.read_csv(file)`` and ``MProVictor.add_category(postera)``.
        :return:
        """
        url = "https://raw.githubusercontent.com/postera-ai/" + \
              "COVID_moonshot_submissions/master/covid_submissions_all_info.csv"
        s = requests.get(url).content
        postera = pd.read_csv(io.StringIO(s.decode('utf-8')))
        cls.add_category(postera)
        return postera

    @classmethod
    def add_category(cls, postera: pd.DataFrame) -> None:
        """
        Postera table has categories as True/False. But it is unlikely that there are multiple.
        Turns out these categories are **not** user submitted.
        However, for consistency with other analysis by other people these are used.

        :param postera: pandas table modified in place
        :return:
        """
        def get_category(row):
            for category in ('Acrylamide', 'Chloroacetamide', 'Vinylsulfonamide', 'Nitrile'):
                if row[category] in ('True', 'true', True):
                    return category
            else:
                return 'non-covalent'

        postera['category'] = postera.apply(get_category, axis=1)

    @classmethod
    def analyse_postera(cls):
        pass

    @classmethod
    def from_postera_row(cls, row: pd.Series, results:Optional=None):
        if row.fragments == 'x0072' or str(row.fragments) == 'nan':
            # these are not hit inspired.
            cls.journal.error(f'No valid inspiration hits for {row.CID}.')
            return None
        elif results and row.CID in results:
            cls.journal.info(f'{row.CID} has already been done.')
            return None
        elif row.covalent_warhead in (False, 'False', 'false'):
            # parse
            return cls.from_hit_codes(hit_codes=row.fragments.split(','),
                                      category='noncolavent')\
                      .place(long_name=row.CID,
                             smiles=row.SMILES)
        elif row.category not in ('Acrylamide', 'Chloroacetamide', 'Vinylsulfonamide', 'Nitrile'):
            cls.journal.warning(f'What is {row["CID"]}? Treating like a non-covalent.')
            return cls.from_hit_codes(hit_codes=row.fragments.split(','),
                                      category='noncolavent')\
                      .place(long_name=row.CID,
                             smiles=row.SMILES)
        else:
            return cls.from_hit_codes(hit_codes=row.fragments.split(','),
                                      category=row.category)\
                      .place(long_name=row.CID,
                             smiles=row.SMILES)

#### Update the dataset with constraints specific for MPro

for cname, con in [('chloroacetamide', 'AtomPair  H  145A  OY  1B HARMONIC 2.1 0.2\n'),
                   ('nitrile', 'AtomPair  H  145A  NX  1B HARMONIC 2.1 0.2\n'),
                   ('acrylamide', 'AtomPair  H  143A  OZ  1B HARMONIC 2.1 0.2\n'),
                   ('vinylsulfonamide', 'AtomPair  H  143A  OZ1 1B HARMONIC 2.1 0.2\n'),
                    ('bromoalkyne', 'AtomPair  H  145A  CY  1B HARMONIC 2.1 0.2\n'),
                   ]:
    MProVictor.add_constraint_to_warhead(name=cname, constraint=con)






