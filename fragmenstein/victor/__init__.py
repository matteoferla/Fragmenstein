from __future__ import annotations

########################################################################################################################

__doc__ = \
    """
??
And can be called via the class method ``.laboratory`` where he can process multiple compounds at once.

The class Victor is assembled via a various class that inherit VictorCommon
This is made by series of classes, whose dependency is not really linear, 
but for simplicity is have been written that way.

1. VictorBase => constructor
2. VictorSafety => error catching
3. VictorJournal => logging
4. VictorPlonk => place
5. VictorOverridables => empty methods
5. VictorStore => save
6. VictorIgor => call igor
7. VictorCommon => common

* VctorPlace => placement
* VictorCombine => merging/linking
* VictorUtils => Bits and bobs
* 
    """
########################################################################################################################

from ._victor_utils import _VictorUtils
from ._victor_validate import _VictorValidate
from ._victor_combine import _VictorCombine
from ._victor_place import _VictorPlace


from .minimalPDB import MinimalPDBParser


class Victor(_VictorUtils, _VictorValidate, _VictorCombine, _VictorPlace):
    """
    Victor (after Dr Victor Frankenstein) is a class that uses both Monster (makes blended compounds)
    and Igor (energy minimises).
    This master reanimator keeps a ``.journal`` (logging, class attribute).

    The constructor sets the protein detail. While, place or combine deal do the analyses.

    :ivar apo_pdbblock: The apo protein template PDB block (inputted)
    :vartype apo_pdbblock: str
    :ivar atomnames: an optional dictionary that gets used by ``Params.from_smiles`` to assigned atom names (inputted)
    :vartype atomnames: Union[None, List, Dict]
    :ivar category: MPro only.
    :vartype category: None/str
    :cvar constrained_atoms: number of atoms constrained (dynamic)
    :vartype constrained_atoms: int
    :ivar constraint: constrains object from rdkit_to_params
    :vartype constraint: Constraints
    :cvar constraint_function_type: name of constraint function. Best not to change.
    :vartype constraint_function_type: str
    :cvar covalent_definitions: definitions of warheads (advanced)
    :vartype covalent_definitions: list
    :ivar covalent_resi: the residue index (PDB) for the covalent attachment (int or int + chain) or reference residue
    :vartype covalent_resi: str
    :ivar covalent_resn:  reference residue name3. the residue name for the covalent attachment.
        For now can only be 'CYS' (or anything else if not covalent)
    :vartype covalent_resn: str
    :ivar energy_score: dict of splits of scores
    :vartype energy_score: dict
    :ivar error_msg: error message if an error of the type error_to_catch was raised and caught
    :vartype error_msg: str
    :cvar error_to_catch: catch error_to_catch.
    :vartype error_to_catch: tuple
    :ivar extra_constraint: extra constraint text
    :vartype extra_constraint: str
    :ivar hits:
    :vartype hits: List[Chem.Mol]
    :ivar igor: igor object
    :vartype igor: Igor
    :ivar is_covalent:
    :vartype is_covalent: bool/None
    :ivar joining_cutoff: max distance between joining mol
    :vartype joining_cutoff: float
    :cvar journal: log
    :vartype journal: Logger
    :ivar ligand_resi:  the residue index (PDB) for the ligand.
    :vartype ligand_resi: str
    :ivar ligand_resn: the residue name for the ligand.
    :vartype ligand_resn: str
    :ivar long_name: name for files
    :vartype long_name: str
    :ivar merging_mode:
    :vartype merging_mode: str
    :ivar minimized_mol:
    :vartype minimized_mol: Mol
    :ivar minimized_pdbblock:
    :vartype minimized_pdbblock: str
    :ivar mmerging_mode:
    :vartype mmerging_mode: str
    :ivar modifications:
    :vartype modifications: dict
    :ivar mol:
    :vartype mol: Mol
    :ivar monster:
    :vartype monster: Monster
    :cvar monster_average_position:
    :vartype monster_average_position: bool
    :cvar monster_mmff_minisation:
    :vartype monster_mmff_minisation: bool
    :cvar monster_throw_on_discard:
    :vartype monster_throw_on_discard: bool
    :ivar mrmsd:
    :vartype mrmsd: mRSMD
    :ivar params:
    :vartype params: Params
    :ivar pose_fx:
    :vartype pose_fx: function
    :cvar possible_definitions:
    :vartype possible_definitions: list
    :cvar quick_reanimation:
    :vartype quick_reanimation: bool
    :ivar reference_mol:
    :vartype reference_mol: NoneType
    :ivar smiles:
    :vartype smiles: str
    :ivar tick:
    :vartype tick: float
    :ivar tock:
    :vartype tock: float
    :ivar unbound_pose:
    :vartype unbound_pose: Pose
    :cvar unconstrained_heavy_atoms:
    :vartype unconstrained_heavy_atoms: int
    :ivar unminimized_pdbblock:
    :vartype unminimized_pdbblock: str
    :cvar warhead_definitions:
    :vartype warhead_definitions: list
    :ivar warhead_harmonisation:
    :vartype warhead_harmonisation: str
    :cvar work_path: class attr. where to save stuff
    :vartype work_path: str

    ``warhead_definitions`` and ``covalent_definitions`` are class attributes that can be modified beforehand to
    allow a new attachment. ``covalent_definitions`` is a list of dictionaries of 'residue', 'smiles', 'names',
    which are needed for the constraint file making. Namely smiles is two atoms and the connection and names is the
    names of each. Cysteine is ``{'residue': 'CYS', 'smiles': '*SC', 'names': ['CONN3', 'SG', 'CB']}``.
    While ``warhead_definitions`` is a list of 'name' (name of warhead for humans),
    'covalent' (the smiles of the warhead, where the zeroth atom is the one attached to the rest),
    'noncovalent' (the warhead unreacted),
    'covalent_atomnames' and 'noncovalent_atomnames' (list of atom names).
    The need for atomnames is actually not for the code but to allow lazy tweaks and analysis downstream
    (say typing in pymol: `show sphere, name CX`).
    """
    pass




