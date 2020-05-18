import logging


class _VictorBaseMixin:
    fragmenstein_merging_mode = 'partial'
    fragmenstein_debug_draw = False
    hits_path = 'hits'
    work_path = 'output'
    journal = logging.getLogger('Fragmenstein')
    journal.setLevel(logging.DEBUG)

    covalent_definitions = [{'residue': 'CYS', 'smiles': '*SC', 'atomnames': ['CONN3', 'SG', 'CB']}]
    warhead_definitions = [{'name': 'nitrile',
                            'covalent': 'C(=N)*',  # zeroth atom is attached to the rest
                            'covalent_atomnames': ['CX', 'NX', 'CONN1'],
                            'noncovalent': 'C(#N)',  # zeroth atom is attached to the rest
                            'noncovalent_atomnames': ['CX', 'NX']
                            },
                           {'name': 'acrylamide',
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
                           {'name': 'vinylsulfonamide',
                            'covalent': 'S(=O)(=O)CC*',  # the N may be secondary etc. so best not do mad substitutions.
                            'covalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX', 'CONN1'],  # OZ tauto
                            'noncovalent': 'S(=O)(=O)C=C',
                            'noncovalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX']
                            }
                           ]

    _connected_names = ('CONN', 'LOWE', 'UPPE', 'CONN1', 'CONN2', 'CONN3', 'LOWER', 'UPPER')

    def __init__(self):
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
        self.fragmenstein = 'Fragmenstein'
        self.unminimised_pdbblock = str()
        self.igor = 'Igor'
        self.minimised_pdbblock = str()
