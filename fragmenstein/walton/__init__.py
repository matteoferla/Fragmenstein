from ._advmoves import WaltonAdvMove
from ._art import WaltonArt
from ._base import WaltonBase
from ._movements import WaltonMove


class Walton(WaltonAdvMove, WaltonArt):
    """
    This is a diagnostic/demo class.
    It aligns molecules for merging to test a hypothesis.
    In the Frankenstein novel,
    (Captain Walton is told the story of the creation of Adam by Victor Frankenstein.
    Basically he is just a spectator but actually allows the story to be known:
    this class does nothing bar demonstrate what can be done)

    Given a list of mols (or smiles as named arguments via the ``from_smiles`` classmethod),
    It aligns them via the ``align_by_mcs`` or ``align_by_map`` methods.
    Calling the instance, will merge molecules and yield the merged molecule.
    The initial molecules are stored as a list in the attribute ``mols``.
    The merged molecule is stored in the attribute ``merged``.
    """

