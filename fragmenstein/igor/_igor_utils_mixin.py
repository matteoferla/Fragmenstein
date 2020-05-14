########################################################################################################################

__doc__ = \
    """
These are extra functionality for Igor
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################

import requests, shutil, pyrosetta
from typing import Optional


class _IgorUtilsMixin:

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
