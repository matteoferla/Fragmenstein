# fragmenstein package

## Subpackages


* fragmenstein.cli package


    * Module contents


* fragmenstein.igor package


    * Module contents


* fragmenstein.laboratory package


    * Submodules


    * fragmenstein.laboratory.laboratory module


    * fragmenstein.laboratory.make_pyrosetta_options module


    * Module contents


* fragmenstein.monster package


    * Submodules


    * fragmenstein.monster.bond_provenance module


    * fragmenstein.monster.positional_mapping module


    * fragmenstein.monster.unmerge_mapper module


    * Module contents


* fragmenstein.mpro package


    * Module contents


* fragmenstein.victor package


    * Submodules


    * fragmenstein.victor.minimalPDB module


    * Module contents


## Submodules

## fragmenstein.m_rmsd module

Combined RMSD


### class fragmenstein.m_rmsd.mRSMD(followup, hits, mappings)
Bases: `object`

RMSD are unaligned and in Å.

> The RMSD has been calculated differently.
> The inbuilt RMSD calculations in RDKit (`Chem.rdMolAlign.GetBestRMS`) align the two molecules,
> this does not align them.
> This deals with the case of multiple hits.
> As a comparision, For euclidean distance the square root of the sum of the differences in each coordinates is taken.
> As a comparision, For a regular RMSD the still-squared distance is averaged before taking the root.
> Here the average is done across all the atom pairs between each hit and the followup.
> Therefore, atoms in followup that derive in the blended molecule by multiple atom are scored multiple times.

> $sqrt{

rac{sum_{i}^{N_{
m{hits}}} (sum_{i}^{n} (q_{i,
m{x}} - h_{i,
m{x}})^2 + (q_{i,
m{y}} - h_{i,
m{y}})^2 + (q_{i,
m{z}} - h_{i,
m{z}})^2 }{ncdot m}}$


#### \__init__(followup, hits, mappings)
This is not meant to be called directly.
mappings is a list of len(hits) containing lists of tuples of atom idx that go from followup to hit

The hit _Name must match that in origin!
currected output of monster.origin_from_mol() or cls.get_origins(to-be-scored-mol, annotated)


* **Parameters**

    
    * **followup** (`Mol`) – the followup compounds


    * **hits** (`Sequence`[`Mol`]) – the fragment hits


    * **mappings** (`List`[`List`[`Tuple`[`int`, `int`]]]) – a complicated affair…



#### calculate_msd(molA, molB, mapping)
A nonroot rmsd.


* **Parameters**

    
    * **molA** – 


    * **molB** – 


    * **mapping** – lists of tuples of atom idx that go from molA to molB



* **Returns**

    nonroot rmsd



#### calculate_rmsd(molA, molB, mapping)

#### classmethod copy_all_possible_origins(annotated, target)
Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
However, the atom order seems to be maintained but I dont trust it. Also dummy atoms are stripped.


* **Parameters**

    
    * **annotated** (`Mol`) – 


    * **target** (`Mol`) – 



* **Return type**

    `Tuple`[`List`[`Mol`], `List`[`List`[`int`]]]



* **Returns**

    a list of mols and a list of orgins (a list too)



#### classmethod copy_origins(annotated, target)
Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
However, the atom order seems to be maintained but I dont trust it. Also dummy atoms are stripped.


* **Parameters**

    
    * **annotated** (`Mol`) – 


    * **target** (`Mol`) – 



* **Returns**

    a list of origins



#### classmethod from_annotated_mols(annotated_followup, hits=None)
Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
This classmethod accepts a followup with has this.


* **Parameters**

    
    * **annotated_followup** (`Mol`) – 


    * **hits** (`Optional`[`Sequence`[`Mol`]]) – 



* **Return type**

    `mRSMD`



* **Returns**

    


#### classmethod from_internal_xyz(annotated_followup)
This is an alternative for when the atoms have _x, _y, _z


* **Parameters**

    **annotated_followup** – 



* **Returns**

    


#### classmethod from_other_annotated_mols(followup, hits, annotated)

* **Return type**

    `mRSMD`



#### classmethod from_unannotated_mols(moved_followup, hits, placed_followup)
Mapping is done by positional overlap between placed_followup and hits
This mapping is the applied to moved_followup.


* **Parameters**

    
    * **moved_followup** (`Mol`) – The mol to be scored


    * **hits** (`Sequence`[`Mol`]) – the hits to score against


    * **placed_followup** (`Mol`) – the mol to determine how to score



* **Return type**

    `mRSMD`



* **Returns**

    


#### classmethod is_origin_annotated(mol)

* **Return type**

    `bool`



#### classmethod is_xyz_annotated(mol)

* **Return type**

    `bool`



#### classmethod migrate_origin(mol, tag='_Origin')
The origin list may be saved as a molecule property rather than an atom -saved as a mol say.


* **Parameters**

    
    * **mol** (`Mol`) – mol to fix


    * **tag** – name of prop



* **Return type**

    `Mol`



* **Returns**

    the same mol



#### classmethod mock()
## Module contents

See GitHub documentation
