# fragmenstein.igor package

## Module contents

Igor energy minises the blended compound using pyrosetta.


### class fragmenstein.igor.Igor(pose, constraint_file, ligand_residue='LIG', key_residues=None)
Bases: `fragmenstein.igor._igor_init._IgorInit`, `fragmenstein.igor._igor_min._IgorMin`, `fragmenstein.igor._igor_utils._IgorUtils`

Regular Igor(..) accepts pyrosetta pose.
`Igor.from_pdbblock(..)` accepts pdb block as str,
while `Igorfrom_pdbfile(..)` accepts filename as str.

`ligand` can be one of many things. default is ‘LIG’. But it can be


* pose index (123)
* PDB index ‘123A’
* a tuple of (PDB resi, PDB chain)
* a residue name in uppercase “LIG”


* a pyrosetta.Vector1 where 1 == the ligand.

If key_residues is None, only the connecting residue is added (if present in the LINK record).
This is overridden if one of many options are given.
If it is a pyrosetta.Vector1 it is assumed that 1 mean select this residue (as a result of a `selector.apply(pose)` operation)
If it is a list or tuple, the elements are interpreted similarly to ligand.


#### residues_in_selector(pose, selector)
This method is just for checking purposes for humans basically.


* **Return type**

    `List`[`str`]



#### residues_in_vector(pose, vector)
This method is just for checking purposes for humans basically.


* **Return type**

    `List`[`str`]
