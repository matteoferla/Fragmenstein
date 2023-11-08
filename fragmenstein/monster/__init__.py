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
from .mcs_mapping import IndexMap, SpecialCompareAtoms, ExtendedFMCSMode, transmute_FindMCS_parameters
from ._ff import MinizationOutcome


# ---  Monster ---------------------------------------------------------------------------------------------------------

from ._communal import _MonsterCommunal
from ._ff import _MonsterFF  # _MonsterFF extends _MonsterUtil
from ._combine import _MonsterCombine  # inherits _MonsterCommunal adds the combine method
from ._place import _MonsterPlace  # inherits _MonsterCommunal adds the combine method


class Monster(_MonsterFF, _MonsterPlace, _MonsterCombine):
    """
    This creates a stitched together monster.
    For initilialisation for either placing or combining, it needs a list of hits (rdkit.Chem.Mol).

    Note, the hits have to be 3D embedded as present in the protein —it would defeat the point otherwise!
    For a helper method to extract them from crystal structures see `Victor.extract_mol`.

    The calculation are done either by place or merge.

    ## Place

    >>> monster.place(mol)

    Given a RDKit molecule and a series of hits it makes a spatially stitched together version
    of the initial molecule based on the hits.
    The reason is to do place the followup compound to the hits as faithfully as possible
    regardless of the screaming forcefields.

    * ``.mol_options`` are the possible equiprobable alternatives.
    * ``.positioned_mol`` is the desired output (rdkit.Chem.Mol object)
    * ``.initial_mol`` is the input (rdkit.Chem.Mol object), this is `None` in a `.combine` call.
    * ``.modifications['scaffold']`` is the combined version of the hits (rdkit.Chem.Mol object).
    * ``.modifications['chimera']`` is the combined version of the hits, but with differing atoms made \
                    to match the followup (rdkit.Chem.Mol object).

    ``.get_positional_mapping``, which works also as a class method,
    creates a dictionary of mol_A atom index to mol_B atom index
    based on distance (cutoff 2&Aring;) and not MCS.

    The code works in two broad steps, first a scaffold is made, which is the combination of the hits (by position).
    Then the followup is placed. It is not embedded with constrained embedding functionality of RDKit as this
    requires the reference molecule to have a valid geometry, which these absolutely do not have this.
    Novel side chains are added by aligning an optimised conformer against the closest 3-4 reference atoms.
    Note that ``.initial_mol`` is not touched. ``.positioned_mol`` may have lost some custom properties,
    but the atom indices are the same.

    If an atom in a Chem.Mol object is provided via ``attachment`` argument and the molecule contains a dummy atom.
    Namely element R in mol file or * in string.

    ## Combine

    >>> monster.combine(keep_all=True, collapse_rings=True, joining_cutoff= 5))

    Combines the hits by merging and linking. ``collapse_rings`` argument results in rings being collapsed
    before merging to avoid oddities.
    The last step within the call is fixing any oddities of impossible chemistry via the call ``rectify``.
    This uses the separate class ``Rectifier`` to fix it.

    ## Attributes

    Common input derived

    :ivar hits:
    :vartype hits: list
    :cvar throw_on_discard: filled by keep_all
    :vartype throw_on_discard: bool

    Common derived

    :var matched: (dynamic) accepted hit names
    :vartype matched: List[str]
    :ivar unmatched: discarded hit names
    :vartype unmatched: List[str]

    :cvar journal: The "journal" is the log of Dr Victor Frankenstein (see Victor for more)
    :vartype journal: Logger
    :ivar modifications: copies of the mols along the way
    :vartype modifications: dict
    :ivar mol_options: equally valid alternatives to self.positioned_mol
    :vartype mol_options: list


    ``place`` specific:

    :ivar positioned_mol:
    :vartype positioned_mol: Mol
    :ivar attachment:
    :vartype attachment: NoneType
    :ivar initial_mol:
    :vartype initial_mol: NoneType
    :ivar average_position:
    :vartype average_position: bool
    :var num_common: (dynamic) number of atoms in common between follow-up and hits
    :vartype num_common: int
    :var percent_common: (dynamic) percentage of atoms of follow-up that are present in the hits
    :vartype percent_common: float

    ``combine`` specific:

    :cvar joining_cutoff: how distant (in Å) is too much?
    :vartype joining_cutoff: int
    :cvar atoms_in_bridge_cutoff: how many bridge atoms can be deleted?
        (0 = preserves norbornane, 1 = preserves adamantane)
    :vartype atoms_in_bridge_cutoff: int

    Class attributes best ignored:

    :cvar closeness_weights: list of functions to penalise closeness (ignore for most applications)
    :vartype closeness_weights: list
    :cvar dummy: The virtual atom where the targets attaches. by default `*`. Best not override.
    :vartype dummy: Mol
    :cvar dummy_symbol: The virtual atom where the targets attaches. by default `*`. Best not override.
    :vartype dummy_symbol: str
    :cvar matching_modes:
    :vartype matching_modes: list
    """
    pass
