########################################################################################################################

__doc__ = \
    """
This is a variant of Victor for MPro that uses data from PostEra
    """

########################################################################################################################


from typing import List, Optional
from warnings import warn

import pandas as pd
import pyrosetta
from rdkit import Chem

from ..demo import MPro as data  # noqa this is basically a package, but actually a class
from ..victor import Victor
from .dataframe import fetch_postera, read_postera


class MProVictor(Victor):
    constraint_function_type = 'FLAT_HARMONIC'

    @classmethod
    def from_hit_codes(cls, hit_codes: List[str], **options):
        hits = [data.get_mol(xnumber) for xnumber in hit_codes]
        return cls(hits=hits, **options)

    @staticmethod
    def pose_fx(pose: pyrosetta.Pose):
        """
        Histidine protonated on epsilon — HIE.
        """
        pdb2pose = pose.pdb_info().pdb2pose
        r = pdb2pose(res=41, chain='A')
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res='HIS').apply(pose)

    @staticmethod
    def poised_pose_fx(pose: pyrosetta.Pose):
        """
        Histidine protonated on delta (HID) and cysteine in thiolate.
        """
        # Americans spell poise as poise not poize.
        pdb2pose = pose.pdb_info().pdb2pose
        r = pdb2pose(res=41, chain='A')
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res='HIS_D').apply(pose)
        r = pdb2pose(res=145, chain='A')
        MutateResidue(target=r, new_res='CYZ').apply(pose)

    def __init__(self, category:Optional[str]=None, **options):
        # this category flag is solely for Mpro?
        # it stems from the Moonshot file.
        self.category = category
        if category == 'noncolavent':
            fx = self.__class__.poised_pose_fx
        else:
            fx = self.__class__.pose_fx
        # --------------------
        defaults = dict(
            pdb_block=data.get_template(),
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

    def place(self, *args, **options):
        defaults = self._determine_extras(options['smiles'] if 'smiles' in options else args[0])
        # singledispatchmethod does not like named args
        kwargs = {**defaults, **options}
        return super().place(*args, **kwargs)

    # ======= postera csv file ops =====================================================================================

    @classmethod
    def fetch_postera(cls):
        __doc__ = fetch_postera.__doc__
        warn('method moved to mpro_data (as function)', category=DeprecationWarning)
        return fetch_postera()

    @classmethod
    def analyse_postera(cls):
        raise NotImplementedError('Incomplete functionality')

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
            cls.journal.info(f'What is {row["CID"]}? Treating like a non-covalent.')
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






