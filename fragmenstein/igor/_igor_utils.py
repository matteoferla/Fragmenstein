########################################################################################################################

__doc__ = \
    """
These are extra functionality for Igor
    """

########################################################################################################################

import requests, shutil
from typing import Optional, Dict, Union
from .pyrosetta_import import pyrosetta  # the real mcCoy or a mock.
from ._igor_base import _IgorBase
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
from enum import Enum


class _IgorUtils(_IgorBase):

    def dock(self) -> pyrosetta.Pose:
        """
        Docks the pose the normal way and without constraints.

        :return:
        """
        docked = self.pose.clone()
        # self.ligand_residue[0] is pose index. This is a bit weird / hacky for now:
        # the pdb_info needs resetting everything chain A.
        docked.pdb_info().set_resinfo(res=self.ligand_residue[0], chain_id='B', pdb_res=1)
        docked.remove_constraints()
        pyrosetta.rosetta.protocols.docking.setup_foldtree(docked, 'A_B', pyrosetta.Vector1([1]))
        scorefxn = pyrosetta.create_score_function('ligand')
        docking = pyrosetta.rosetta.protocols.docking.DockMCMProtocol()
        docking.set_scorefxn(scorefxn)
        docking.apply(docked)
        return docked

    def per_atom_scores(self,
                        pose: Optional[pyrosetta.Pose]=None,
                        target_res: Optional[int]=None,
                        scorefxn: Optional[pyrosetta.ScoreFunction]=None) -> Dict[str, Dict[str, float]]:
        """
        Per atom scores are generally a bad idea as a score relative to something else is better.
        So do treat with the appropriate caution.
        NB. these scores will not sum to the that of the residue.

        Given a pose, a target_res and a scorefxn return a dict of per atom scores:

        * Lenard-Jones attraction 6-term
        * Lenard-Jones repulsion 12-term
        * Solvatation (zero)
        * Electrostatic interactions

        :param pose:
        :param target_res:
        :param scorefxn:
        :return: a dict of atom names to dict of 'lj_atr', 'lj_rep', 'fa_solv', 'fa_elec' to value
        """
        # Defaults
        if target_res is None:
            target_res = self.ligand_residue[0]
        if pose is None:
            pose = self.pose
        if scorefxn is None:
            # The minimisation is done under cartesian setting though.
            # and constrained!
            scorefxn = pyrosetta.get_fa_scorefxn()
        # Prep
        score_types = ['lj_atr', 'lj_rep', 'fa_solv', 'fa_elec']
        residue = pose.residue(target_res)
        scores = {residue.atom_name(i): {st: 0 for st in score_types} for i in range(1, residue.natoms() + 1)}
        # Iterate per target residue's atom per all other residues' atoms
        for i in range(1, residue.natoms() + 1):
            iname = residue.atom_name(i)
            for r in range(1, pose.total_residue() + 1):
                other = pose.residue(r)
                for o in range(1, other.natoms() + 1):
                    score = pyrosetta.toolbox.atom_pair_energy.etable_atom_pair_energies(residue,
                                                                                         i,
                                                                                         other,
                                                                                         o,
                                                                                         scorefxn)
                    for st, s in zip(score_types, score):
                        # technically only fa_intra_rep is 100 fold less than fa_rep
                        scores[iname][st] += s if target_res != r else s/100
        return scores

    def display_energy(self,
                       minimized_mol: Union[Chem.Mol, None],
                       term: Union[str, Enum]='all',
                       scorefxn: Optional[pyrosetta.ScoreFunction]=None,
                       colorMap:str='coolwarm') -> Dict[str, float]:
        """
        Given the ``minimized_mol`` display the energy contributions.
        As Victor has this while Igor is naive of bond orders, if None is passed it will made due with no bond order.

        ``term`` is the term to use:

        * ``all`` will be ``['lj_atr', 'lj_rep', 'fa_solv', 'fa_elec']``
        * ``lenard-jones`` will ``['lj_atr', 'lj_rep']``
        * while ``['lj_atr', 'lj_rep', 'fa_solv', 'fa_elec']`` will be the individual ones.

        As this is rather obscure, there is an enum as an attribute of ``display_energy`` method called ``Term``,
        with the values ``.display.Term.ALL``, ``.display.Term.ATTRACTION``, ``.display.Term.REPULSION``,
        ``.display.Term.LENARDJONES``, ``.display.Term.SOLVATATION``, ``.display.Term.ELECTROSTATIC``
        which is possibly more sane (e.g. ``igor.display_energy(minimized_mol, term=igor.display.Term.ELECTROSTATIC)``).
        """
        if minimized_mol is None:
            minimized_mol:Chem.Mol = self.mol_from_pose(add_dummy=False)
        if scorefxn is None:
            scorefxn:pyrosetta.ScoreFunction = pyrosetta.get_fa_scorefxn()
        # ['all', 'lenard-jones', 'lj_atr', 'lj_rep', 'fa_solv', 'fa_elec']
        # deal with choices of terms:
        if isinstance(term, Enum):
            term = term.value
        weights = {'lj_atr': scorefxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_atr),
                   'lj_rep': scorefxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_rep),
                   'fa_solv': scorefxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_sol),
                   'fa_elec': scorefxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_elec)}
        if term == 'all':
            pass
        elif term == 'lenard-jones':
            weights['fa_solv'] = 0
            weights['fa_elec'] = 0
        elif term in weights:
            for other in set(weights.keys()) - {term}:
                weights[other] = 0
        else:
            raise ValueError(f'Term {term} is not `all`,`lj_atr`,`lj_rep`,`fa_solv`,`fa_elec`')

        # get scores
        scores = self.per_atom_scores()

        def get_combined(atom: Chem.Atom):
            # combine the scores
            scored = scores[atom.GetPDBResidueInfo().GetName()]
            return sum([scored[k] * weights[k] for k in scored])

        contribs = [get_combined(atom) for atom in minimized_mol.GetAtoms()]
        # make into a 2D mol
        mol = Chem.Mol(minimized_mol)
        AllChem.Compute2DCoords(mol)
        # displays automatically:
        SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap=colorMap, contourLines=10)
        return {atom.GetProp('atom_name'): round(get_combined(atom), 1) for atom in minimized_mol.GetAtoms()}

    @classmethod
    def download_map(cls, pdbcode: str, filename: str):
        """
        Download to disk a CCP4 map of a given PDB code
        :param pdbcode:
        :param filename:
        :return:
        """
        assert '.ccp4' in filename, f'This downloads ccp4 maps ({filename})'
        r = requests.get(f'http://www.ebi.ac.uk/pdbe/coordinates/files/{pdbcode.lower()}.ccp4', stream=True)
        if r.status_code == 200:
            with open(filename, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)

    @classmethod
    def relax_with_ED(cls, pose,
                      ccp4_file: str,
                      constraint_file: Optional[str] = None,
                      constraint_weights=(30, 20, 10),
                      cycles: int = 5
                      ) -> pyrosetta.Pose:
        """
        Relaxes ``pose`` based on the ccp4 electron density map provided. See ``download_map`` to download one.

        :param pose:
        :param ccp4_file: download map from ePDB
        :param constraint_weights: weights for the constraints First run with high weights, then lower.
        :param cycles: number of cycles
        :return: Relaxes pose in place
        """
        scorefxnED = pyrosetta.get_fa_scorefxn()
        ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(ccp4_file)
        sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
        sdsm.apply(pose)
        ## Set ED constraint
        elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
        scorefxnED.set_weight(elec_dens_fast, constraint_weights[0])
        ## Relax
        for w in (30, 20, 10):
            cls.fix_scorefxn_weights(scorefxn=scorefxnED,
                                     elec_dens_fast=w,
                                     atom_pair_constraint=w,
                                     angle_constraint=w,
                                     dihedral_constraint=w)
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, cycles)
            relax.apply(pose)
        return pose

    @classmethod
    def add_constraint_file(cls, pose: pyrosetta.Pose, constraint_file: str):
        """
        Add a constraint file to the pose.

        :param pose:
        :param constraint_file:
        :return:
        """
        stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
        setup.constraint_file(constraint_file)
        setup.apply(pose)

    @classmethod
    def fix_scorefxn_weights(cls, scorefxn, **constraint_weights):
        for contype_name, weight in constraint_weights.items():
            contype = stm.score_type_from_name(contype_name)
            scorefxn.set_weight(contype, weight)

    @classmethod
    def tethered_relax(cls,
                       original: pyrosetta.Pose,
                      constraint_weight: float=5,
                      cycles: int=15,
                      relax_to_start_coords:bool=True) -> pyrosetta.Pose:
        """
        Relax by constraining the original pose.

        :param original: the pose to relax, but makes a copy
        :param constraint_weight: the weight of the constraints
                                (`coordinate_constraint`, `angle_constraint`, `atom_pair_constraint`)
        :param cycles: number of cycles
        :param relax_to_start_coords: relax to the start coords or not
        :return: the relaxed pose

        ``relax_to_start_coords`` is the same as ``constraint_weight=0`` if the original pose has no custom constraints
        """
        pose: pyrosetta.Pose = original.clone()
        # configure constraints
        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
        scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, constraint_weight)
        scorefxn.set_weight(pr_scoring.ScoreType.angle_constraint, constraint_weight)
        scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, constraint_weight)
        pyrosetta.rosetta.basic.options.set_boolean_option('relax:constrain_relax_to_start_coords', relax_to_start_coords)
        pyrosetta.rosetta.basic.options.set_boolean_option('relax:coord_constrain_sidechains', relax_to_start_coords)
        # set up the relax sampler
        pyrosetta.rosetta.protocols.relax.FastRelax.register_options()
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.constrain_relax_to_start_coords(relax_to_start_coords)
        relax.apply(pose)
        return pose

    @staticmethod
    def init_pyrosetta():
        """
        Common options for pyrosetta. Does not use `pyrosetta_help` atm
        """
        pyrosetta.init(extra_options='-no_optH false -ex1 -ex2 -mute all -ignore_unrecognized_res false ' +\
                                     '-load_PDB_components false -ignore_waters true')


class Term(Enum):
    ALL = 'all'
    ATTRACTION = 'lj_atr'
    REPULSION = 'lj_rep'
    LENARDJONES = 'lenard-jones'
    SOLVATATION = 'fa_solv'
    ELECTROSTATIC = 'fa_elec'

setattr(_IgorUtils.display_energy, 'Term', Term)