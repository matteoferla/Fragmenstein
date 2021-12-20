########################################################################################################################

__doc__ = \
    """
Base metho=ods
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################


import logging, warnings
import re
import time
from typing import List, Union, Optional, Callable, Dict

from rdkit import Chem
from ..m_rmsd import mRSMD
from ..monster import Monster


class _VictorBase:
    quick_reanimation = False  # thorugh reanimation?
    monster_average_position = False
    monster_throw_on_discard = False
    monster_mmff_minisation = True
    constraint_function_type = 'FLAT_HARMONIC'
    work_path = 'output'
    journal = logging.getLogger('Fragmenstein')
    journal.setLevel(logging.DEBUG)

    covalent_definitions = [{'residue': 'CYS', 'smiles': '*SC', 'atomnames': ['CONN3', 'SG', 'CB']}]
    warhead_definitions = [{'name': 'acrylamide',
                            'covalent': 'C(=O)CC*',  # the N may be secondary etc. so best not do mad substitutions.
                            'covalent_atomnames': ['CZ', 'OZ', 'CY', 'CX', 'CONN1'],
                            # OZ needs to tautomerise & h-bond happily.
                            'noncovalent': 'C(=O)C=C',
                            'noncovalent_atomnames': ['CZ', 'OZ', 'CY', 'CX']},
                           {'name': 'chloroacetamide',
                            'covalent': 'C(=O)C*',  # the N may be secondary etc. so best not do mad substitutions.
                            'covalent_atomnames': ['CY', 'OY', 'CX', 'CONN1'],
                            # OY needs to tautomerise & h-bond happily.
                            'noncovalent': 'C(=O)C[Cl]',
                            'noncovalent_atomnames': ['CY', 'OY', 'CX', 'CLX']
                            },
                           {'name': 'nitrile',
                            'covalent': 'C(=N)*',  # zeroth atom is attached to the rest
                            'covalent_atomnames': ['CX', 'NX', 'CONN1'],
                            'noncovalent': 'C(#N)',  # zeroth atom is attached to the rest
                            'noncovalent_atomnames': ['CX', 'NX']
                            },
                           {'name': 'vinylsulfonamide',
                            'covalent': 'S(=O)(=O)CC*',  # the N may be secondary etc. so best not do mad substitutions.
                            'covalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX', 'CONN1'],  # OZ tauto
                            'noncovalent': 'S(=O)(=O)C=C',
                            'noncovalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX']
                            },
                           {'name': 'bromoalkyne',
                            'covalent': 'C(=C)*',
                            'covalent_atomnames': ['CX', 'CY', 'CONN1'],
                            # OY needs to tautomerise & h-bond happily.
                            'noncovalent': 'C#C[Br]',
                            'noncovalent_atomnames': ['CX', 'CY', 'BRX']
                            },
                           ]

    # these may be wrong and need checking.
    possible_definitions = [{'name': 'aurothiol',  # gold salt
                             'covalent': 'S[Au]*',
                             'covalent_atomnames': ['SY', 'AUX', 'CONN1'],
                             # OY needs to tautomerise & h-bond happily.
                             'noncovalent': 'S[Au]P(CC)(CC)CC',
                             'noncovalent_atomnames': ['SY', 'AUX', 'PL', 'CL1', 'CL2', 'CL3', 'CL4', 'CL5', 'CL6']
                             },
                            {'name': 'aldehyde',
                             'covalent': 'C(O)*',
                             'covalent_atomnames': ['CX', 'OX', 'CONN1'],
                             'noncovalent': '[C:H1]=O',  # this at
                             'noncovalent_atomnames': ['CX', 'OX']
                             },
                            ]

    _connected_names = ('CONN', 'LOWE', 'UPPE', 'CONN1', 'CONN2', 'CONN3', 'LOWER', 'UPPER')

    error_to_catch = ()  # Exception

    remove_other_hetatms = True #remove all protein heteroatoms that are not water or ligand when loaded from PDB

    # ================== Init ==========================================================================================

    def __init__(self,
                 hits: List[Chem.Mol],
                 pdb_filename: Union[None, str] = None,
                 pdb_block: Union[None, str] = None,
                 ligand_resn: str = 'LIG',
                 ligand_resi: Union[int, str] = '1B',
                 covalent_resn: str = 'CYS',  # no other option is accepted.
                 covalent_resi: Optional[Union[int, str]] = None,
                 extra_protein_constraint: Union[str] = None,
                 pose_fx: Optional[Callable] = None,
#                 remove_other_hetatms: bool = True,
                 ):
        """
        Initialise Victor in order to allow either combinations (merging/linking without a given aimed for molecule)
        or placements (using a given aimed for molecule).

        :param hits: list of rdkit molecules
        :param pdb_filename: file of apo structure
        :param pdb_block: alternative for above: a string of apo structure
        :param ligand_resn: 3 letter code or your choice
        :param ligand_resi: Rosetta-style pose(int) or pdb(str)
        :param covalent_resn: only CYS accepted. if smiles has no * it is ignored
        :param covalent_resi: Rosetta-style pose(int) or pdb(str)
        :param extra_protein_constraint: multiline string of constraints relevant to the protein
        :param pose_fx: a function to call with pose to tweak or change something before minimising.

        """
        # ## Store
        # entry attributes
        if pdb_filename:
            with open(pdb_filename) as fh:
                self.apo_pdbblock = fh.read()
        elif pdb_block:
            self.apo_pdbblock = pdb_block
        else:
            raise ValueError('Provide a pdb_filename or pdb_block of the template')
        self.hits = hits
        self.ligand_resn = ligand_resn.upper()
        self.ligand_resi = ligand_resi
        self.covalent_resn = covalent_resn.upper()
        self.covalent_resi = covalent_resi
        self._correct_covalent_resi()  # defined in plonk. todo: split into covalent and anchor residue.
        self.extra_constraint = extra_protein_constraint
        self.pose_fx = pose_fx
        # ## Fill by place and combine differently
        self.long_name = 'ligand'
        self.smiles = None
        # ## Filled by place
        self.merging_mode = "none_permissive"
        # ## Filled by combine
        self.joining_cutoff = None
        # ## Calculated
        self.is_covalent = None
        self.params = None
        self.mol = None
        self.constraint = None
        self.modifications = {}
        self.unminimized_pdbblock = None
        self.monster = Monster(hits,
                               average_position=self.monster_average_position)
        self.igor = None
        self.unbound_pose = None
        self.minimized_pdbblock = None
        self.minimized_mol = None
        self.reference_mol = None  # filled only for validate
        # buffers etc.
        self._warned = []
        self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                             'unbound_ref2015': {'total_score': float('nan')}}
        self.mrmsd = mRSMD.mock()
        self.ddG = float('nan')
        # for debug purposes
        self.tick = time.time()
        self.tock = float('inf')
        self.error_msg = ''

    # ----------------- init called methods ----------------------------------------------------------------------------
    @classmethod
    def slugify(cls, name: str):
        return re.sub(r'[\W_.-]+', '-', name)
