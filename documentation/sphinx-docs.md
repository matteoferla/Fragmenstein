# fragmenstein.core package

## Module contents

This is Fragmenstein proper. and contains the class `Fragmenstein`.


### class fragmenstein.core.Fragmenstein(mol, hits, attachment=None, debug_draw=False)
Bases: `fragmenstein.core._utility_mixin._FragmensteinUtil`

Given a RDKit molecule and a series of hits it makes a spatially stitched together version of the initial molecule based on the hits.
The reason is to do place the followup compound to the hits as faithfully as possible regardless of the screaming forcefields.


* `.scaffold` is the combined version of the hits (rdkit.Chem.Mol object).


* `.chimera` is the combined version of the hits, but with differing atoms made to match the followup (rdkit.Chem.Mol object).


* `.positioned_mol` is the desired output (rdkit.Chem.Mol object)

Note, the hits have to be spatially aligned —i.e. extracted from crystal structures in bond form.

`.get_positional_mapping`, which works also as a class method, creates a dictionary of mol_A atom index to mol_B atom index
based on distance (cutoff 2&Aring;) and not MCS.

The code works in two broad steps, first a scaffold is made, which is the combination of the hits (by position).
Then the followup is placed. It is not embedded with constraint embed, which requires the reference molecule to have a valid geometry.
`.scaffold` and `.chimera` and `.positioned_mol` absolutely do not have this.
Novel side chains are added by aligning an optimised conformer against the closest 3-4 reference atoms.
Note that `.initial_mol` is not touched. `.positioned_mol` may have lost some custom properties, but the atom idices are the same.

If an atom in a Chem.Mol object is provided via `attachment` argument and the molecule contains a dummy atom as
defined in the `dummy` class variable. Namely element R in mol file or \* in string is the default.


#### \__init__(mol, hits, attachment=None, debug_draw=False)
Initialize self.  See help(type(self)) for accurate signature.


#### dummy( = <rdkit.Chem.rdchem.Mol object>)
The virtual atom where the targets attaches


#### dummy_symbol( = '\*')

#### get_mcs_mapping(molA, molB)
This is a weird method. It does a strict MCS match.
And then it uses laxer searches and finds the case where a lax search includes the strict search.


* **Parameters**

    
    * **molA** – query molecule


    * **molB** – target/ref molecule



* **Return type**

    `Tuple`[`Dict`[`int`, `int`], `dict`]



* **Returns**

    mapping and mode



#### classmethod get_positional_mapping(mol_A, mol_B, cutoff=2)
Returns a map to convert overlapping atom of A onto B
Cutoff 2 &Aring;.


* **Parameters**

    
    * **mol_A** (`Mol`) – first molecule (Chem.Mol) will form keys


    * **mol_B** (`Mol`) – second molecule (Chem.Mol) will form values



* **Return type**

    `Dict`[`int`, `int`]



* **Returns**

    dictionary mol A atom idx -> mol B atom idx.



#### make_chimera()
This is to avoid extreme corner corner cases. E.g. here the MCS is ringMatchesRingOnly=True and AtomCompare.CompareAny,
while for the positioning this is not the case.


* **Return type**

    `Mol`



* **Returns**

    


#### merge(scaffold, fragmentanda, anchor_index, attachment_details)

* **Return type**

    `Mol`



#### merge_hits()
Recursively stick the hits together and average the positions.


* **Return type**

    `Mol`



* **Returns**

    the rdkit.Chem.Mol object that will fill `.scaffold`



#### place_from_map(mol=None)

* **Return type**

    `Mol`



#### posthoc_refine(scaffold)
Averages the overlapping atoms.


* **Parameters**

    **scaffold** – 



* **Returns**

    


#### pretweak()
What if the fragments were prealigned slightly? Really bad things.


* **Return type**

    `None`



* **Returns**

# fragmenstein.igor package

## Module contents

Igor energy minises the blended compound using pyrosetta.


### class fragmenstein.igor.Igor(pose, constraint_file, ligand_residue='LIG', key_residues=None)
Bases: `fragmenstein.igor._igor_init_mixin._IgorInitMixin`, `fragmenstein.igor._igor_min_mixin._IgorMinMixin`, `fragmenstein.igor._igor_utils_mixin._IgorUtilsMixin`

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

# fragmenstein.victor package

## Module contents

Victor (after Dr Victor Frankenstein) is a class that uses both Fragmenstein (makes blended compounds) and Igor (energy minimises).
This master reanimator keeps a `.journal` (logging, class attribute).
And can be called via the class method `.laboratory` where he can process multiple compounds at once.


### class fragmenstein.victor.Victor(smiles, hits, pdb_filename, long_name='ligand', ligand_resn='LIG', ligand_resi='1B', covalent_resn='CYS', covalent_resi=None, extra_constraint=None, pose_fx=None)
Bases: `object`


* `smiles` SMILES string (inputted)


* `long_name` name for files


* `ligand_resn` the residue name for the ligand.


* `ligand_resi` the residue index (PDB) for the ligand.


* `covalent_resi` the residue index (PDB) for the covalent attachment


* `covalent_resn` the residue name for the covalent attachment. For now can only be ‘CYS’


* `params` Params instance


* `constraint` Constraint or None depending on if covalent.


* `mol` the molecule


* `covalent_definitions` class attr. that stores for each possible attachment residue (CYS) defs for constraints.


* `warhead_definitions` class attr. that stores warheader info


* `journal` class attr. logging


* `work_path` class attr. where to save stuff

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
Adding a ‘constraint’ to an entry will apply that constraint.


#### \__init__(smiles, hits, pdb_filename, long_name='ligand', ligand_resn='LIG', ligand_resi='1B', covalent_resn='CYS', covalent_resi=None, extra_constraint=None, pose_fx=None)

* **Parameters**

    
    * **smiles** (`str`) – smiles of followup, optionally covalent (_e.g._ `\*CC(=O)CCC`)


    * **hits** (`List`[`Mol`]) – list of rdkit molecules


    * **pdb_filename** (`str`) – file of apo structure


    * **long_name** (`str`) – gets used for filenames so will get slugified


    * **ligand_resn** (`str`) – 3 letter code or your choice


    * **ligand_resi** (`Union`[`int`, `str`]) – Rosetta-style pose(int) or pdb(str)


    * **covalent_resn** (`str`) – only CYS accepted. if smiles has no \* it is ignored


    * **covalent_resi** (`Union`[`int`, `str`, `None`]) – Rosetta-style pose(int) or pdb(str)


    * **extra_constraint** (`Optional`[`str`]) – multiline string of constraints..


    * **pose_fx** (`Optional`[`Callable`]) – a function to call with pose to tweak or change something before minimising.



#### classmethod add_constraint_to_warhead(name, constraint)
Add a constraint (multiline is fine) to a warhead definition.
This will be added and run by Igor’s minimiser.


* **Parameters**

    
    * **name** (`str`) – 


    * **constraint** (`str`) – 



* **Returns**

    None



#### classmethod closest_hit(pdb_filenames, target_resi, target_chain, target_atomname, ligand_resn='LIG')
This classmethod helps choose which pdb based on which is closer to a given atom.


* **Parameters**

    
    * **pdb_filenames** (`List`[`str`]) – 


    * **target_resi** (`int`) – 


    * **target_chain** (`str`) – 


    * **target_atomname** (`str`) – 


    * **ligand_resn** – 



* **Return type**

    `str`



* **Returns**

    


#### classmethod copy_names(acceptor_mol, donor_mol)

#### covalent_definitions( = [{'residue': 'CYS', 'smiles': '\*SC', 'atomnames': ['CONN3', 'SG', 'CB']}])

#### classmethod enable_logfile(filename='reanimation.log', level=20)
The journal is output to a file.


* **Parameters**

    
    * **filename** – file to write.


    * **level** – logging level



* **Return type**

    `None`



* **Returns**

    None



#### classmethod enable_stdout(level=20)
The `cls.journal` is output to the terminal.


* **Parameters**

    **level** – logging level



* **Return type**

    `None`



* **Returns**

    None



#### hits_path( = 'hits')

#### journal( = <Logger Fragmenstein (DEBUG)>)

#### classmethod laboratory(entries, cores=1)

#### classmethod make_all_warhead_combinations(smiles, warhead_name)
Convert a unreacted warhead to a reacted one in the SMILES


* **Parameters**

    
    * **smiles** (`str`) – unreacted SMILES


    * **warhead_name** (`str`) – name in the definitions



* **Return type**

    `Optional`[`str`]



* **Returns**

    SMILES



#### classmethod make_covalent(smiles, warhead_name=None)
Convert a unreacted warhead to a reacted one in the SMILES


* **Parameters**

    
    * **smiles** (`str`) – unreacted SMILES


    * **warhead_name** (`Optional`[`str`]) – name in the definitions. If unspecified it will try and guess (less preferrable)



* **Return type**

    `Optional`[`str`]



* **Returns**

    SMILES



#### pose_mod_step()
This method is intended for make inherited mods easier.
:return:


#### post_igor_step()
This method is intended for make inherited mods easier.
:return:


#### post_fragmenstein_step()
This method is intended for make inherited mods easier.
:return:


#### post_params_step()
This method is intended for make inherited mods easier.
:return:


#### classmethod slack_me(msg)
Send message to a slack webhook


* **Parameters**

    **msg** (`str`) – Can be dirty and unicode-y.



* **Returns**

    did it work?



* **Return type**

    bool



#### slugify(name)

#### warhead_definitions( = [{'name': 'nitrile', 'covalent': 'C(=N)\*', 'covalent_atomnames': ['CX', 'NX', 'CONN1'], 'noncovalent': 'C(#N)', 'noncovalent_atomnames': ['CX', 'NX']}, {'name': 'acrylamide', 'covalent': 'C(=O)CC\*', 'covalent_atomnames': ['CZ', 'OZ', 'CY', 'CX', 'CONN1'], 'noncovalent': 'C(=O)C=C', 'noncovalent_atomnames': ['CZ', 'OZ', 'CY', 'CX']}, {'name': 'chloroacetamide', 'covalent': 'C(=O)C\*', 'covalent_atomnames': ['CY', 'OY', 'CX', 'CONN1'], 'noncovalent': 'C(=O)C[Cl]', 'noncovalent_atomnames': ['CY', 'OY', 'CX', 'CLX']}, {'name': 'vinylsulfonamide', 'covalent': 'S(=O)(=O)CC\*', 'covalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX', 'CONN1'], 'noncovalent': 'S(=O)(=O)C=C', 'noncovalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX']}])

#### work_path( = 'output')

## fragmenstein.m_rmsd module


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
currected output of fragmenstein.origin_from_mol() or cls.get_origins(to-be-scored-mol, annotated)


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

#### classmethod copy_origins(annotated, target)
Fragmenstein leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
However, the atom order seems to be maintained but I dont trust it. Also dummy atoms are stripped.


* **Parameters**

    
    * **annotated** (`Mol`) – 


    * **target** (`Mol`) – 



* **Returns**

    a list of origins



#### classmethod from_annotated_mols(annotated_followup, hits)
Fragmenstein leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
This classmethod accepts a followup with has this.


* **Parameters**

    
    * **annotated_followup** (`Mol`) – 


    * **hits** (`Sequence`[`Mol`]) – 



* **Returns**

    


#### classmethod from_other_annotated_mols(followup, hits, annotated)

#### classmethod from_unannotated_mols(moved_followup, hits, placed_followup)