# fragmenstein.victor package

## Submodules

## fragmenstein.victor.minimalPDB module


### class fragmenstein.victor.minimalPDB.MinimalPDBParser(block)
Bases: `object`

This purpose build PDB parser simply fixes the serial numbers.
The reason is that writing a custom 50 line class is easier that
having biopython or other non-builtin requirement as a requirement
Importing the PDB into RDKit is inadvisable.


#### \__init__(block)
Initialize self.  See help(type(self)) for accurate signature.


#### append(other)
Add a second parser data to it. But only its coordinates and connections.


#### get_max_serial()

* **Return type**

    `int`



#### get_serial(entry)

#### offset_connections(offset)

* **Return type**

    `None`



#### offset_serials(offset)

* **Return type**

    `None`



#### parse(block)

* **Return type**

    `None`



#### set_serial(entry, value)

* **Return type**

    `None`


## Module contents


The class Victor is assembled via a various class that inherit VictorCommon
This is made by series of classes, whose dependency is not really linear, 
but for simplicity is have been written that way.

1. VictorBase => constructor
2. VictorSafety => error catching
3. VictorJournal => logging
4. VictorPlonk => place
5. VictorOverridables => empty methods
5. VictorStore => save
6. VictorIgor => call igor
7. VictorCommon => common

* VctorPlace => placement
* VictorCombine => merging/linking
* VictorUtils => Bits and bobs


### class fragmenstein.victor.Victor(hits, pdb_filename, ligand_resn='LIG', ligand_resi='1B', covalent_resn='CYS', covalent_resi=None, extra_protein_constraint=None, pose_fx=None)
Bases: `fragmenstein.victor._victor_utils._VictorUtils`, `fragmenstein.victor._victor_validate._VictorValidate`, `fragmenstein.victor._victor_combine._VictorCombine`, `fragmenstein.victor._victor_place._VictorPlace`

and Igor (energy minimises).
This master reanimator keeps a `.journal` (logging, class attribute).

The constructor sets the protein detail. While, place or combine deal do the analyses.


* **Variables**

    
* **apo_pdbblock** (*str*) – The apo protein template PDB block (inputted)
* **atomnames** (*Union**[**None**, **List**, **Dict**]*) – an optional dictionary that gets used by `Params.from_smiles` to assigned atom names (inputted)
* **category** (*None/str*) – MPro only.
* **constrained_atoms** (*int*) – number of atoms constrained (dynamic)
* **constraint** (*Constraints*) – constrains object from rdkit_to_params
* **constraint_function_type** (*str*) – name of constraint function. Best not to change.
* **covalent_definitions** (*list*) – definitions of warheads (advanced)
* **covalent_resi** (*str*) – the residue index (PDB) for the covalent attachment (int or int + chain) or reference residue
* **covalent_resn** (*str*) – reference residue name3. the residue name for the covalent attachment.
    
For now can only be ‘CYS’ (or anything else if not covalent)

* **energy_score** (*dict*) – dict of splits of scores
* **error_msg** (*str*) – error message if an error of the type error_to_catch was raised and caught
* **error_to_catch** (*tuple*) – catch error_to_catch.
* **extra_constraint** (*str*) – extra constraint text
* **hits** (*List**[**Chem.Mol**]*) – 
* **igor** (*Igor*) – igor object
* **is_covalent** (*bool/None*) – 
* **joining_cutoff** (*float*) – max distance between joining mol
* **journal** (*Logger*) – log
* **ligand_resi** (*str*) – the residue index (PDB) for the ligand.
* **ligand_resn** (*str*) – the residue name for the ligand.
* **long_name** (*str*) – name for files
* **merging_mode** (*str*) – 
* **minimised_mol** (*Mol*) – 
* **minimised_pdbblock** (*str*) – 
* **modifications** (*dict*) – 
* **mol** (*Mol*) – 
* **monster** (*Monster*) – 
* **monster_average_position** (*bool*) – 
* **monster_mmff_minisation** (*bool*) – 
* **monster_throw_on_discard** (*bool*) – 
* **mrmsd** (*mRSMD*) – 
* **params** (*Params*) – 
* **pose_fx** (*function*) – 
* **possible_definitions** (*list*) – 
* **quick_renanimation** (*bool*) – 
* **reference_mol** (*NoneType*) – 
* **smiles** (*str*) – 
* **tick** (*float*) – 
* **tock** (*float*) – 
* **unbound_pose** (*Pose*) – 
* **unconstrained_heavy_atoms** (*int*) – 
* **unminimised_pdbblock** (*str*) – 
* **warhead_definitions** (*list*) – 
* **warhead_harmonisation** (*str*) – 
* **work_path** (*str*) – class attr. where to save stuff

`warhead_definitions` and `covalent_definitions` are class attributes that can be modified beforehand to
allow a new attachment. `covalent_definitions` is a list of dictionaries of ‘residue’, ‘smiles’, ‘names’,
which are needed for the constraint file making. Namely smiles is two atoms and the connection and names is the
names of each. Cysteine is `{'residue': 'CYS', 'smiles': '\*SC', 'names': ['CONN3', 'SG', 'CB']}`.
While `warhead_definitions` is a list of ‘name’ (name of warhead for humans),
‘covalent’ (the smiles of the warhead, where the zeroth atom is the one attached to the rest),
‘noncovalent’ (the warhead unreacted),
‘covalent_atomnames’ and ‘noncovalent_atomnames’ (list of atom names).
The need for atomnames is actually not for the code but to allow lazy tweaks and analysis downstream
(say typing in pymol: show sphere, name CX).
