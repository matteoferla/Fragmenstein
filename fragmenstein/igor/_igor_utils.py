########################################################################################################################

__doc__ = \
    """
These are extra functionality for Igor
    """

########################################################################################################################

import requests, shutil, pyrosetta
from typing import Optional, Dict

import pyrosetta


class _IgorUtils:

    def dock(self) -> pyrosetta.Pose:
        """
        Docks the pose the normal way and without constraints.

        :return:
        """
        docked = self.pose.clone()
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
                        scores[iname][st] += s
        return scores

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
    def relax_with_ED(cls, pose, ccp4_file: str, constraint_file: Optional[str] = None) -> None:
        """
        Relaxes ``pose`` based on the ccp4 electron density map provided. See ``download_map`` to download one.
        :param pose:
        :param ccp4_file: download map from ePDB
        :return: Relaxes pose in place
        """
        scorefxnED = pyrosetta.get_fa_scorefxn()
        ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(ccp4_file)
        sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
        sdsm.apply(pose)
        ## Set ED constraint
        elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
        scorefxnED.set_weight(elec_dens_fast, 30)
        ## Set generic constraints
        if constraint_file:
            stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
            for contype_name in ("atom_pair_constraint", "angle_constraint", "dihedral_constraint"):
                contype = stm.score_type_from_name(contype_name)
                scorefxnED.set_weight(contype, 5)
            setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
            setup.constraint_file(constraint_file)
            setup.apply(pose)
        ## Relax
        for w in (30, 20, 10):
            scorefxnED.set_weight(elec_dens_fast, w)
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 5)
            relax.apply(pose)
