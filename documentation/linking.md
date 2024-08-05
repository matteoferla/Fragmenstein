## Linking

The linker atoms are typically a carbon linear chain, where the first atom is nitrogen or oxygen in versions prior to 1.1, or set by the user.

## Method
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

## Cutoff

The cutoff can be changed in one of four ways:

* the environment variant `FRAGMENSTEIN_CUTOFF` ,
* `settings['cutoff']`,
* the argument in Monster initialisation `joining_cutoff` or 
* the Victor class attribute `monster_joining_cutoff`

### Linking

The word 'joining' is used in the code as it did not initial link with new atoms, so was in the grey area between linking and merging.
The linking is done by the private method `_join_atoms`, which will add atoms in a line if needed.
Currently, the added linker is a hydrocarbon.
The distance between two atoms added is 1.22 as when minimised they will be in a zip-zap.


### Linker atom

In version prior to 1.1 the first linker atom was oxygen, this was changed to N.

The identity is controlled by the class attribute `linker_atom` in Monster,
or the Victor argument `linker_atom` (or via env variable `FRAMENSTEIN_LINKER_ATOM`).
Element-specific atomic radii are not considered, so changing to tellurium will behave as if it were oxygen.

 The special rule for the first atom is to better model a likely non-hydrophobic element at that position. Initially, this was oxygen but was changed to nitrogen because ether is not a common reaction moiety, whereas a bridging nitrogen is a common product of several reactions, such as reductive amination. 
 
The default linker length is short because catalogue/reaction-product enumeration is preferred for longer linkers due to patchiness in catalogue space and issues with excessive rotatable bonds (as rigidification incurs an entropic penalty). 

Fragmenstein is nevertheless helpful in such cases. 
For example, I, Matteo, perform a search with a dot-separated degenerate SMARTS pattern derived from the parent hits and places the compounds with Fragmenstein (https://github.com/matteoferla/Arthorian-Quest).
In Wills et al. 2023 (https://pubs.acs.org/doi/10.1021/acs.jcim.3c00276), a search for superstructures that include two parent hits was conducted and then placed with Fragmenstein.

## Caveat

Parent fragment hits that are disconnected will be connected.
If your crystal structure is a homodimer, make sure the hit used is a single copy from the correct chain.