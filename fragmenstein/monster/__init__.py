########################################################################################################################
__doc__ = \
    """
This is Monster proper. and contains the class ``Monster``.
The inheritance is as follows:


**_MonsterCombine**  inherits:

* ``_MonsterMerge``  (``_MonsterCommunal`` <- ``_MonsterTracker`` <- ``_MonsterBase``)
* ``_MonsterRing``

**_MonsterPlace**  inherits:

* ``_MonsterMerge``  (``_MonsterBlend`` < - ``_MonsterCommunal`` <- ``_MonsterTracker`` <- ``_MonsterBase``)

**_MonsterUtils** inherits

* ``_MonsterCommunal``  ( <- ``_MonsterTracker`` <- ``_MonsterBase``) 
* ``GPM``

Where:

**_MonsterBase** adds the cvars and the ``__init__`` and its dependent methods

**_MonsterTracker** inherits ``_MonsterBase`` and adds just a method to better store modifications

**_MonsterCommunal** inherits ``_MonsterTracker`` ( <- ``_MonsterBase``)


It adds methods for misc purposes beyond place/combine.
It is inherited by ``Monster`` only.
    """

########################################################################################################################

# --- independent of Monster -------------------------------------------------------------------------------------------
from .positional_mapping import GPM
from .bond_provenance import BondProvenance
from .unmerge_mapper import Unmerge


# ---  Monster ---------------------------------------------------------------------------------------------------------

from ._communal import _MonsterCommunal
from ._utility import _MonsterUtil  # Adds extras, not called by place/combine
from ._combine import _MonsterCombine  # inherits _MonsterCommunal adds the combine method
from ._place import _MonsterPlace  # inherits _MonsterCommunal adds the combine method


class Monster(_MonsterUtil, _MonsterPlace, _MonsterCombine):
    """
    This creates a stitched together monster.
    For initilialisation for either placing or merging, it needs a list of hits (rdkit.Chem.Mol).

    The calculation are done either by place or merge.

    ## Place
    monster.place(mol)
    Given a RDKit molecule and a series of hits it makes a spatially stitched together version of the initial molecule based on the hits.
    The reason is to do place the followup compound to the hits as faithfully as possible regardless of the screaming forcefields.

    * ``.scaffold`` is the combined version of the hits (rdkit.Chem.Mol object).
    * ``.mol_options`` are the possible mol_options to use.
    * ``.chimera`` is the combined version of the hits, but with differing atoms made to match the followup (rdkit.Chem.Mol object).
    * ``.positioned_mol`` is the desired output (rdkit.Chem.Mol object)

    Note, the hits have to be spatially aligned —i.e. extracted from crystal structures in bond form (see. `extract_mol`).

    ``.get_positional_mapping``, which works also as a class method, creates a dictionary of mol_A atom index to mol_B atom index
    based on distance (cutoff 2&Aring;) and not MCS.

    The code works in two broad steps, first a scaffold is made, which is the combination of the hits (by position).
    Then the followup is placed. It is not embedded with constraint embed, which requires the reference molecule to have a valid geometry.
    ``.scaffold`` and ``.chimera`` and ``.positioned_mol`` absolutely do not have this.
    Novel side chains are added by aligning an optimised conformer against the closest 3-4 reference atoms.
    Note that ``.initial_mol`` is not touched. ``.positioned_mol`` may have lost some custom properties, but the atom idices are the same.

    If an atom in a Chem.Mol object is provided via ``attachment`` argument and the molecule contains a dummy atom as
    defined in the ``dummy`` class variable. Namely element R in mol file or * in string is the default.
    """
    pass
