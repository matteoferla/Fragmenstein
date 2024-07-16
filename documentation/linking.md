## Linking

The linking of two unconnected compounds happens in Monster class's `_join_neighboring.py` methods.

The method `join_neighboring_mols` does the joining.

    monster.join_neighboring_mols(mol_A, mol_B) # -> mol_combined

### Choosing

The closest atoms are chosen by the private method `_find_closest`.
This method depends on the `closeness_weights` class-attribute, a list of tuples of `function(Chem.Atom) -> bool` and float penalty.
Namely, if the passed atom is True the penalty is added to the distance.

* warhead atom: NaN (don't touch them)
* fully bonded: +1.0
* ring atom: +0.5

### Joining

The joining is done by the private method `_join_atoms`, which will add atoms in a line.
The distance between two atoms added is 1.22 as when minimised they will be in a zip-zap.
Currently, the added linker is a hydrocarbon.

### Linker atom

In version prior to 1.1 the first linker atom was oxygen, this was changed to N.

The identity is controlled by the class attribute `linker_atom` in Monster,
or the Victor argument `linker_atom` (or via env variable `FRAMENSTEIN_LINKER_ATOM`).
Element-specific atomic radii are not considered, so changing to tellurium will behave as if it were oxygen.