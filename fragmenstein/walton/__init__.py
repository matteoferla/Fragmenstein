from ._advmoves import WaltonAdvMove
from ._art import WaltonArt
from ._base import WaltonBase
from ._movements import WaltonMove
from ._polygon import WaltonPolygon


class Walton(WaltonAdvMove, WaltonArt, WaltonPolygon):
    """
    This is a diagnostic/demo class.
    It superposes molecules for merging to test a hypothesis.
    In the Frankenstein novel,
    (Captain Walton is told the story of the creation of Adam by Victor Frankenstein.
    Basically he is just a spectator but actually allows the story to be known:
    this class does nothing bar demonstrate what can be done)

    Given a list of mols (or smiles as named arguments via the ``from_smiles`` classmethod),
    It superposes them via the ``superpose_by_mcs`` or ``superpose_by_map`` methods.
    Calling the instance, will merge molecules and yield the merged molecule.
    The initial molecules are stored as a list in the attribute ``mols``.
    The merged molecule is stored in the attribute ``merged``.

    Walton has some methods to speed up operations

    .. code-block:: python

        demo = Walton.from_smiles(spiropentadiene='C1=CC21C=C2')
        demo.ring_on_plane(mol_idx=0, ring_idx=1)
        demo.print_coords(mol_idx=0)

    Walton is aimed at making life easier to set up an artificial merger condition and show it.

    .. code-block:: python

        demo = Walton.from_smiles(resorcinol='c1ccc(O)cc1O', eugenol='Oc1ccc(cc1OC)CC=C') # create instance
        demo.superpose_by_map({(0,1):{4:0, 3:1, 2:2}})  # superpose molecules by atom indices
        demo()  # merge (Fragmenstein's Monster)
        demo.show3d()  # show 2d, names and 3d

    The latter is equivalent to:

    .. code-block:: python

        display(demo)
        demo.to_nglview()

    The key feature of the class is that it allows easy transformations

    .. code-block:: python

    demo = Walton.from_smiles(furan='c1ccco1')
    display(demo)  # find atom numbers w/o counting
    demo.ring_on_plane(mol_idx=0)  # flat on xy
    demo.atom_on_axis(mol_idx=0, atom_idx=4, axis='x') # ox on axis
    # the atom_idx can be actually be passed a Geometry.Point3D:
    demo.atom_on_axis(mol_idx=0,
                     atom_idx=demo.get_centroid_of_atoms(1, 2, mol_idx=0),
                    axis='x')
    """
