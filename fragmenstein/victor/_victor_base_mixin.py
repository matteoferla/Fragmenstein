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


import logging


class _VictorBaseMixin:
    quick_renanimation = False  # thorugh reanimation?
    adam_merging_mode = 'none_permissive'
    adam_debug_draw = False
    adam_average_position = False
    adam_joining_cutoff = 5.  # Å
    adam_throw_on_discard = False
    adam_mmff_minisation = True
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
    possible_definitions =[{'name': 'aurothiol', # gold salt
                            'covalent': 'S[Au]*',
                            'covalent_atomnames': ['SY', 'AUX', 'CONN1'],
                            # OY needs to tautomerise & h-bond happily.
                            'noncovalent': 'S[Au]P(CC)(CC)CC',
                            'noncovalent_atomnames': ['SY', 'AUX', 'PL', 'CL1', 'CL2', 'CL3', 'CL4', 'CL5', 'CL6']
                            },
                           {'name': 'aldehyde',
                            'covalent': 'C(O)*',
                            'covalent_atomnames': ['CX', 'OX', 'CONN1'],
                            'noncovalent': '[C:H1]=O', # this at
                            'noncovalent_atomnames': ['CX', 'OX']
                            },
                           ]

    _connected_names = ('CONN', 'LOWE', 'UPPE', 'CONN1', 'CONN2', 'CONN3', 'LOWER', 'UPPER')

    error_to_catch = Exception

    def __init__(self):
        # raise NotImplementedError('Abstract method')
        # gets overridden
        self.long_name = str()
        self.smiles = str()
        self.apo_pdbblock = str()
        self.hits = list()
        self.ligand_resn = str()
        self.ligand_resi = str()
        self.covalent_resn = str()
        self.covalent_resi = str()
        self.extra_constraint = str()
        self.pose_fx = lambda x: None
        self.is_covalent = False
        # these are calculated.
        self.params = 'Params'
        self.mol = 'Chem.Mol'
        self.constraint = 'Constraint'
        self.fragmenstein = 'Adam'
        self.modifications = []
        self.unminimised_pdbblock = str()
        self.igor = 'Igor'
        self.minimised_pdbblock = str()
