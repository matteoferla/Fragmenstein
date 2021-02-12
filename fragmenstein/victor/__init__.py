from __future__ import annotations

########################################################################################################################

__doc__ = \
    """
Victor (after Dr Victor Frankenstein) is a class that uses both Monster (makes blended compounds) and Igor (energy minimises).
This master reanimator keeps a ``.journal`` (logging, class attribute).
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
    * ``smiles`` SMILES string (inputted)
    * ``long_name`` name for files
    * ``ligand_resn`` the residue name for the ligand.
    * ``ligand_resi`` the residue index (PDB) for the ligand.
    * ``covalent_resi`` the residue index (PDB) for the covalent attachment
    * ``covalent_resn`` the residue name for the covalent attachment. For now can only be 'CYS'
    * ``params`` Params instance
    * ``constraint`` Constraint or None depending on if covalent.
    * ``mol`` the molecule
    * ``covalent_definitions`` class attr. that stores for each possible attachment residue (CYS) defs for constraints.
    * ``warhead_definitions`` class attr. that stores warheader info
    * ``journal`` class attr. logging
    * ``work_path`` class attr. where to save stuff

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
    Adding a 'constraint' to an entry will apply that constraint.
    ``merging_mode:str`` is class attributes that control Monster.

    """
    pass




