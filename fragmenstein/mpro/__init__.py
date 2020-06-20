########################################################################################################################

__doc__ = \
    """
This is a variant of Victor for MPro that uses data from PostEra
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################


from ..victor import Victor

import os, pyrosetta
from rdkit import Chem
from typing import List, Optional

import pandas as pd
import io
import requests


def pose_fx(pose):
    """
    Histidine in delta.
    """
    pose2pdb = pose.pdb_info().pdb2pose
    r = pose2pdb(res=41, chain='A')
    MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
    MutateResidue(target=r, new_res='HIS').apply(pose)


def poised_pose_fx(pose):
    """
    Histidine in delta and cysteine in thiolate.
    """
    pose2pdb = pose.pdb_info().pdb2pose
    r = pose2pdb(res=41, chain='A')
    MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
    MutateResidue(target=r, new_res='HIS_D').apply(pose)
    r = pose2pdb(res=145, chain='A')
    MutateResidue(target=r, new_res='CYZ').apply(pose)

class MProVictor(Victor):
    fragmenstein_merging_mode = 'none_permissive'
    constraint_function_type = 'FLAT_HARMONIC'

    @classmethod
    def get_mpro_path(cls):
        if os.path.islink(__file__):
            path = os.readlink(__file__)
        else:
            path = __file__
        return os.path.join(os.path.dirname(path), 'data')

    @classmethod
    def get_mol(cls, xnumber):
        mpro_folder = cls.get_mpro_path()
        xnumber = xnumber.strip()
        mol = Chem.MolFromMolFile(os.path.join(mpro_folder, 'hit_mols', f'Mpro-{xnumber}.mol'))
        mol.SetProp('_Name', xnumber)
        return mol

    @classmethod
    def from_hit_codes(cls, smiles: str, hit_codes:List[str], long_name:str, category:Optional[str]=None):
        hits = [cls.get_mol(xnumber) for xnumber in hit_codes]
        return cls(smiles=smiles, hits=hits, long_name=long_name, category=category)

    def __init__(self, smiles: str, hits:List[Chem.Mol], long_name:str, category:Optional[str]=None):
        mpro_folder = self.get_mpro_path()
        apo = os.path.join(mpro_folder, 'template.pdb')
        atomnames = {}
        if category == 'noncolavent':
            fx = poised_pose_fx
        else:
            fx = pose_fx
        extra_constraint = 'AtomPair  SG  145A  NE2  41A HARMONIC 3.5 0.2\n'
        if category not in (None, 'noncovalent') and '_' in category:
            cname, rxd = category.split('_')
            if rxd == 'noncovalent':
                wd = [wd for wd in self.warhead_definitions if wd['name'] == cname][0]
                mol = Chem.MolFromSmiles(smiles)
                nc = Chem.MolFromSmiles(wd['noncovalent'])
                atomnames = dict(zip(mol.GetSubstructMatch(nc), wd['noncovalent_atomnames']))
                fx = poised_pose_fx
                extra_constraint += 'AtomPair  SG  145A  CX   1B HARMONIC 3.2 0.5\n'
                extra_constraint += wd['constraint']

        super().__init__(smiles=smiles,
                         hits=hits,
                         pdb_filename=apo,
                         long_name=long_name,
                         ligand_resn='LIG',
                         ligand_resi='1B',
                         covalent_resn='CYS', covalent_resi='145A',
                         extra_constraint=extra_constraint,
                         pose_fx=fx,
                         atomnames=atomnames)

    @classmethod
    def combine_codes(cls, hit_codes: List[str], warhead_harmonisation='first'):
        hits = [cls.get_mol(xnumber) for xnumber in hit_codes]
        return cls.combine(hits=hits, warhead_harmonisation=warhead_harmonisation)


    @classmethod
    def combine(cls, hits:List[Chem.Mol], warhead_harmonisation='first'):
        mpro_folder = cls.get_mpro_path()
        apo = os.path.join(mpro_folder, 'template.pdb')
        atomnames = {}
        fx = pose_fx
        extra_constraint = 'AtomPair  SG  145A  NE2  41A HARMONIC 3.5 0.2\n'
        return super().combine(hits=hits,
                         pdb_filename=apo,
                         ligand_resn='LIG',
                         ligand_resi='1B',
                         covalent_resn='CYS', covalent_resi='145A',
                         extra_constraint=extra_constraint,
                         pose_fx=fx,
                         atomnames=atomnames,
                         warhead_harmonisation=warhead_harmonisation)

    # ======= postera csv file ops =====================================================================================

    @classmethod
    def fetch_postera(cls):
        """
        Reads the submission file off Github.
        For a local version, just ``postera = pd.read_csv(file)`` and ``MProVictor.add_category(postera)``.
        :return:
        """
        url = "https://raw.githubusercontent.com/postera-ai/COVID_moonshot_submissions/master/covid_submissions_all_info.csv"
        s = requests.get(url).content
        postera = pd.read_csv(io.StringIO(s.decode('utf-8')))
        cls.add_category(postera)
        return postera

    @classmethod
    def add_category(cls, postera):
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

#### Update the dataset with constraints specific for MPro

for cname, con in [('chloroacetamide', 'AtomPair  H  145A  OY  1B HARMONIC 2.1 0.2\n'),
                   ('nitrile', 'AtomPair  H  145A  NX  1B HARMONIC 2.1 0.2\n'),
                   ('acrylamide', 'AtomPair  H  143A  OZ  1B HARMONIC 2.1 0.2\n'),
                   ('vinylsulfonamide', 'AtomPair  H  143A  OZ1 1B HARMONIC 2.1 0.2\n'),
                    ('bromoalkyne', 'AtomPair  H  145A  CY  1B HARMONIC 2.1 0.2\n'),
                   ]:
    MProVictor.add_constraint_to_warhead(name=cname, constraint=con)






