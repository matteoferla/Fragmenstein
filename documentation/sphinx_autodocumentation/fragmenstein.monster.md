# fragmenstein.monster package

## Submodules

## fragmenstein.monster.bond_provenance module


### class fragmenstein.monster.bond_provenance.BondProvenance()
Bases: `enum.Enum`

Where does the bond come from. This is used to keep names consistent…
For now original is used in places. The others are interchangeable TBH.

* ORIGINAL( = 1)
* MAIN_NOVEL( = 2)
* OTHER_NOVEL( = 3)
* LINKER( = 4)
* UNASSIGNED( = 5)

#### classmethod copy_bond(donor, acceptor)

* **Return type**

    `None`



#### classmethod get_bond(bond)

* **Return type**

    `BondProvenance`



#### classmethod get_bonds(bonds)

* **Return type**

    `List`[`BondProvenance`]



#### classmethod has_bond(bond)

* **Return type**

    `bool`



#### classmethod set_all_bonds(mol, provenance_name)
Sets the provenance of all bonds in mol to a category, which is a string from the provenance


* **Parameters**

    * **mol** (`Mol`) – 
    * **provenance_name** (`str`) – A string original | main_novel ” other_novel | linker



* **Return type**

    `None`



* **Returns**

    


#### classmethod set_bond(bond, provenance_name)

* **Return type**

    `None`



#### classmethod set_bonds(bonds, provenance_name)

* **Return type**

    `None`


## fragmenstein.monster.positional_mapping module

Positional mapping


### class fragmenstein.monster.positional_mapping.GPM()
Bases: `object`

This class simply contains `get_positional_mapping` and is inherited both by Monster and Unmerge.
`get_positional_mapping` teturns a map to convert overlapping atom of A onto B


#### cutoff( = 2)

#### classmethod get_positional_mapping(mol_A, mol_B, dummy_w_dummy=True)
Returns a map to convert overlapping atom of A onto B
Cutoff 2 &Aring; (see class attr.)


**Parameters**
  
* **mol_A** (`Mol`) – first molecule (Chem.Mol) will form keys
* **mol_B** (`Mol`) – second molecule (Chem.Mol) will form values
* **dummy_w_dummy** – match 



* **Return type**

    `Dict`[`int`, `int`]



* **Returns**

    dictionary mol A atom idx -> mol B atom idx.


## fragmenstein.monster.unmerge_mapper module

Unmerge mapper (not inherited)


### class fragmenstein.monster.unmerge_mapper.Unmerge(followup, mols, maps, no_discard=False)
Bases: `fragmenstein.monster.positional_mapping.GPM`

This class tries to solve the mapping problem by try all possible mappings of the target to the ligand.
It is one of three in Monster (full merge, partial merge, unmerge.

It is great with fragments that do not connect, but is bad when a hit has a typo.


* the positions must overlap if any atom is mapped in two maps


* no bond can be over 3 A

The chosen map `combined_map` is a dict that goes from `followup` mol to `combined` mol which
is the hits in a single molecule.

Note that some molecules are discarded entirely.


#### \__init__(followup, mols, maps, no_discard=False)

**Parameters**

    
* **followup** (*Chem.Mol*) – the molecule to place
* **mols** (*List**[**Chem.Mol**]*) – 3D molecules
* **maps** (*Dict**[**List**[**Dict**[**int**, **int**]**]**]*) – can be generated outseide of Monster by `.make_maps`.
* **no_discard** (`bool`) – do not allow any to be discarded

#### bond(idx=None)

#### check_possible_distances(other, possible_map, combined, combined_map, cutoff=3)

#### distance_cutoff( = 3)
how distance is too distant in Å


#### get_inter_distance(molA, molB, idxA, idxB)

* **Return type**

    `float`



#### get_key(d, v)
Given a value and a dict and a value get the key.
:type d: `dict`
:param d:
:type v: `Any`
:param v:
:return:


#### get_possible_map(other, label, o_map, inter_map, combined, combined_map)
This analyses a single map (o_map) and returns a possible map


* **Parameters**

    
* **other** (`Mol`) – 
* **label** (`str`) – 
* **o_map** (`Dict`[`int`, `int`]) – followup -> other
* **inter_map** (`Dict`[`int`, `int`]) – 
* **combined** (`Mol`) – 
* **combined_map** (`Dict`[`int`, `int`]) – followup -> combined


* **Return type**

    `Dict`[`int`, `int`]



* **Returns**

    followup -> other



#### judge_n_move_on(combined, combined_map, other, possible_map, others, disregarded)
The mutables need to be within their own scope


* **Parameters**

    
* **combined** – 
* **combined_map** – 
* **other** – 
* **possible_map** – 
* **others** – 
* **disregarded** – 



* **Returns**

    


#### classmethod make_maps(target, mols, mode=None)
This is basically if someone is using this class outside of Monster

Returns a dictionary of key mol name and
value a list of possible dictionary with idex of target to the index given mol.
Note that a bunch of mapping modes can be found in Monster init mixin class.


* **Parameters**

    
* **target** (`Mol`) – the molecule to be mapped
* **mols** (`List`[`Mol`]) – the list of molecules with positional data to be mapped to
* **mode** (`Optional`[`Dict`[`str`, `Any`]]) – dict of setting for MCS step



* **Return type**

    `Dict`[`str`, `List`[`Dict`[`int`, `int`]]]



* **Returns**

    


#### max_strikes( = 3)
number of discrepancies tollerated.


#### measure_map(mol, mapping)

* **Parameters**

    
* **mol** (`Mol`) – 
* **mapping** (`Dict`[`int`, `int`]) – followup to comined



* **Return type**

    `array`



* **Returns**

    


#### offness(mol, mapping)
How many bonds are too long?
:type mol: `Mol`
:param mol:
:type mapping: `Dict`[`int`, `int`]
:param mapping:
:rtype: `float`
:return:


#### pick( = 0)

#### rotational_approach( = True)

#### store(combined, combined_map, disregarded)

#### template_sorter_factory(accounted_for)
returns the number of atoms that have not already been accounted for.


* **Return type**

    `Callable`



#### unmerge_inner(combined, combined_map, others, disregarded)
Assesses a combination of maps
rejections: unmapped (nothing maps) / unnovel (adds nothing)


* **Parameters**

    
* **combined** (`Mol`) – 
* **combined_map** (`Dict`[`int`, `int`]) – 
* **others** (`List`[`Mol`]) – 
* **disregarded** (`List`[`Mol`]) – 



* **Return type**

    `None`



* **Returns**

    

## Module contents

This is Monster proper. and contains the class `Monster`.
The inheritance is as follows:

**_MonsterCombine**  inherits:


* `_MonsterMerge`  (`_MonsterCommunal` <- `_MonsterTracker` <- `_MonsterBase`)
* `_MonsterRing`

**_MonsterPlace**  inherits:


* `_MonsterMerge`  (`_MonsterBlend` < - `_MonsterCommunal` <- `_MonsterTracker` <- `_MonsterBase`)

**_MonsterUtils** inherits


* `_MonsterCommunal`  ( <- `_MonsterTracker` <- `_MonsterBase`)
* `GPM`

Where:

**_MonsterBase** adds the cvars and the `__init__` and its dependent methods

**_MonsterTracker** inherits `_MonsterBase` and adds just a method to better store modifications

**_MonsterCommunal** inherits `_MonsterTracker` ( <- `_MonsterBase`)

It adds methods for misc purposes beyond place/combine.
It is inherited by `Monster` only.


### class fragmenstein.monster.Monster(hits, average_position=False)
Bases: `fragmenstein.monster._utility._MonsterUtil`, `fragmenstein.monster._place._MonsterPlace`, `fragmenstein.monster._combine._MonsterCombine`

This creates a stitched together monster.
For initilialisation for either placing or combining, it needs a list of hits (rdkit.Chem.Mol).

Note, the hits have to be 3D embedded as present in the protein —it would defeat the point otherwise!
For a helper method to extract them from crystal structures see Victor.extract_mol.

The calculation are done either by place or merge.

## Place

    monster.place(mol)

Given a RDKit molecule and a series of hits it makes a spatially stitched together version
of the initial molecule based on the hits.
The reason is to do place the followup compound to the hits as faithfully as possible
regardless of the screaming forcefields.


* `.mol_options` are the possible equiprobable alternatives.
* `.positioned_mol` is the desired output (rdkit.Chem.Mol object)
* `.initial_mol` is the input (rdkit.Chem.Mol object), this is None in a .combine call.
* `.modifications['scaffold']` is the combined version of the hits (rdkit.Chem.Mol object).
* `.modifications['chimera']` is the combined version of the hits, but with differing atoms made to match the followup (rdkit.Chem.Mol object).

`.get_positional_mapping`, which works also as a class method,
creates a dictionary of mol_A atom index to mol_B atom index
based on distance (cutoff 2&Aring;) and not MCS.

The code works in two broad steps, first a scaffold is made, which is the combination of the hits (by position).
Then the followup is placed. It is not embedded with constrained embedding functionality of RDKit as this
requires the reference molecule to have a valid geometry, which these absolutely do not have this.
Novel side chains are added by aligning an optimised conformer against the closest 3-4 reference atoms.
Note that `.initial_mol` is not touched. `.positioned_mol` may have lost some custom properties,
but the atom indices are the same.

If an atom in a Chem.Mol object is provided via `attachment` argument and the molecule contains a dummy atom.
Namely element R in mol file or \* in string.

## Combine

    monster.combine(keep_all=True, collapse_rings=True, joining_cutoff= 5))

Combines the hits by merging and linking. `collapse_rings` argument results in rings being collapsed
before merging to avoid oddities.
The last step within the call is fixing any oddities of impossible chemistry via the call `rectify`.
This uses the separate class `Rectifier` to fix it.

## Attributes

Common input derived


* **Variables**

    
    * **hits** (*list*) – 
    * **throw_on_discard** (*bool*) – filled by keep_all


Common derived

**Variables**

    
* **matched** (*List**[**str**]*) – (dynamic) accepted hit names
* **unmatched** (*List**[**str**]*) – discarded hit names
* **journal** (*Logger*) – The “journal” is the log of Dr Victor Frankenstein (see Victor for more)
* **modifications** (*dict*) – copies of the mols along the way
* **mol_options** (*list*) – equally valid alternatives to self.positioned_mol

`place` specific:


* **Variables**

    
* **positioned_mol** (*Mol*) – 
* **attachment** (*NoneType*) – 
* **initial_mol** (*NoneType*) – 
* **average_position** (*bool*) – 
* **num_common** (*int*) – (dynamic) number of atoms in common between follow-up and hits
* **percent_common** (*float*) – (dynamic) percentage of atoms of follow-up that are present in the hits


`combine` specific:


* **Variables**

    
* **joining_cutoff** (*int*) – how distant (in Å) is too much?
* **atoms_in_bridge_cutoff** (*int*) – how many bridge atoms can be deleted?
(0 = preserves norbornane, 1 = preserves adamantane)


Class attributes best ignored:


* **Variables**
    
* **closeness_weights** (*list*) – list of functions to penalise closeness (ignore for most applications)
* **dummy** (*Mol*) – The virtual atom where the targets attaches. by default \*. Best not override.
* **dummy_symbol** (*str*) – The virtual atom where the targets attaches. by default \*. Best not override.
* **matching_modes** (*list*) –