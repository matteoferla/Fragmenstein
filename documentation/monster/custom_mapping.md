## Custom mapping for Monster

In the placement operation without prior merging ('none' or 'expand'),
a custom mapping can be supplied to the placement operation
in the form `{'hit1': {1: 2, 3: 4}, 'hit2': {1: 6, 2: 8}}`, where for some hits, 
specified by name (_i.e._ `mol.Prop('_Name')`), is assigned the required correspondence between
hit atom index and followup atom index.

In `monster/mcs_mapping` and `monster/_place_modes/_expand.py` is the code that allows this.
In particular the `SpecialCompareAtoms` class is used to control the MCS comparison.

### SpecialCompareAtoms
The mapping as discussed in GitHub issue #23 is in the format

    mapping = { 'hit1': {1:1,2:5} 'hit2': {3:3,4:4,4:6}}

The hit index is first, followup index is the second.
The index `-1` for a followup index is the same as not providing the hit index,
it is described here solely for clarity not for use.

    mapping = { 'hit1': {1:1,2:5, 3:-1} 'hit2': {3:3,4:4,4:6}}

The index `-2` for a followup index will result in the hit atom index
not matching any followup index.

    mapping = { 'hit1': {1:1,2:5, 3:-2} 'hit2': {3:3,4:4,4:6}}

If ``exclusive_mapping`` argument of `SpecialCompareAtoms.__init__` is True,
then if a followup index is present in one hit, but not in a second hit,
then no atom of the second hit will match that followup atom.
A negative index for a hit atom index means that no atom in that hit will match the
corresponding followup index.

    mapping = { 'hit1': {1:1,2:5,-1:3, -2: 7} 'hit2': {3:3,4:4,4:6}}

However, a positive integer on a different hit overrides it, therefore,
in the above followup atom 3 cannot be matched to any atom in hit1, but will match
atom 3 in hit2. Followup atom 7 will not match either.

### Filters

However, this does not prevent the result being a mapping without those indices.
For that the method `Monster.get_mcs_mappings` makes sure that the mapping does contain those atoms.

### Example

A nice test case comes from an [Angewandte Chemie paper](https://onlinelibrary.wiley.com/cms/asset/fb95854c-3379-4533-912f-d98736bb627c/anie202204052-toc-0001-m.png)
with the following merger:

![merger](https://onlinelibrary.wiley.com/cms/asset/fb95854c-3379-4533-912f-d98736bb627c/anie202204052-toc-0001-m.png)

This is a nice case because it features a troublesome mapping, because the primary hit, F04, 
which looks like a Hantzsch thiazole synthesis product, presents can be flipped either way
upon investigation of the crystal map
and the followup bold4 is indeed a flipped merger. 
Therefore, without a custom_map the followup is incorrect.

Where I to make a flipped thiazole F04 I could do:
```python

with Chem.SDMolSupplier('test_mols/5SB7_mergers.sdf') as reader:
    mols = list(reader)

monster = Monster([mols[0]])
monster.place(Chem.Mol(mols[0]),    # it technically saves a copy already, but one should always be careful
              merging_mode='expansion',
              custom_map={'F04': {-1: 4,  # no amine
                                  4: -2,  # no amine
                                  12: 12,
                                  6: 13,
                                  13: 6,
                                  15: 14,
                                  14: 15}}
              )
monster.show_comparison()
monster.to_nglview()
```
This requires identifying the atom indices.
To find out the atom indices I could do:
```python
monster.draw_nicely(monster.positioned_mol)
```

**Nota bene maxime**: atom indices change when going from a Chem.Mol to SMILES.
In SMILES and mol files the atom indices are as they are read. But when you do the conversion
`Chem.MolToSmiles(mol)` the indices are changed. Therefore, for Victor,
which before May 2022 (v. 0.9) accepted only SMILES, 
one should get the followup index atoms as they appear in the SMILES.
When a followup in the format `Chem.Mol` is passed to `Victor.place` (`funtools.singledispatchmethod` overloaded)
the staticmethod
`Monster.renumber_followup_custom_map` is called to fix these allowing one to provide the indices
as they are in the Chem.Mol, which are renumbered behind the scenes. However, 
calling `Monster.convert_origins_to_custom_map` will return indices in the renumbered format, so cannot be used
as a reference (unless `Monster.renumber_followup_custom_map` is used on it to convert this
from monster.position_mol to the actual original mol).

To find out how the two hits look in 3D, I could do:
```python
monster.to_nglview()
```

```python
from fragmenstein import Igor

Igor.init_pyrosetta()

victor = Victor(hits=mols[:2], pdb_filename='template.pdb', ligand_resi='1X')
              
victor.place(mols[3],
             long_name=mols[3].GetProp('_Name'),
              merging_mode='expansion',
             custom_map={'F36': {1: 7},
                          'F04': {0: -1,  # no amine
                                  -1: 7,  # no amine
                                  13: 5,  # to N
                                  12: 13,  # root to Ph
                                  6: 14,  # to C
                                  }
                          }
             
              )
victor.show_comparison()
victor.to_nglview()
victor.validate(mols[3])['reference2minimized_rmsd'] < 1
```