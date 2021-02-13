# fragmenstein.mpro package

## Module contents

This is a variant of Victor for MPro that uses data from PostEra


### class fragmenstein.mpro.MProVictor(category=None, \*\*options)
Bases: `fragmenstein.victor.Victor`


#### \__init__(category=None, \*\*options)
Initialise Victor in order to allow either combinations (merging/linking without a given aimed for molecule)
or placements (using a given aimed for molecule).


* **Parameters**

    
    * **hits** – list of rdkit molecules


    * **pdb_filename** – file of apo structure


    * **ligand_resn** – 3 letter code or your choice


    * **ligand_resi** – Rosetta-style pose(int) or pdb(str)


    * **covalent_resn** – only CYS accepted. if smiles has no \* it is ignored


    * **covalent_resi** – Rosetta-style pose(int) or pdb(str)


    * **extra_protein_constraint** – multiline string of constraints relevant to the protein


    * **pose_fx** – a function to call with pose to tweak or change something before minimising.



#### classmethod add_category(postera)
Postera table has categories as True/False. But it is unlikely that there are multiple.
Turns out these categories are **not** user submitted.
However, for consistency with other analysis by other people these are used.


* **Parameters**

    **postera** (`DataFrame`) – pandas table modified in place



* **Return type**

    `None`



* **Returns**

    


#### classmethod analyse_postera()

#### constraint_function_type( = 'FLAT_HARMONIC')

#### classmethod fetch_postera()
Reads the submission file off Github.
For a local version, just `postera = pd.read_csv(file)` and `MProVictor.add_category(postera)`.
:return:


#### classmethod from_hit_codes(hit_codes, \*\*options)

#### classmethod from_postera_row(row, results=None)

#### classmethod get_mol(xnumber)

#### classmethod get_mpro_path()

#### place(\*\*options)
Places a followup (smiles) into the protein based upon the hits.
Do note that while Monster’s place accepts a mol, while place_smiles a smiles
Victor’s place accepts only smiles.


* **Parameters**

    
    * **smiles** – smiles of followup, optionally covalent (_e.g._ `\*CC(=O)CCC`)


    * **long_name** – gets used for filenames so will get corrected


    * **merging_mode** – 


    * **atomnames** – an optional dictionary that gets used by `Params.from_smiles`


    * **extra_ligand_constraint** – 



* **Returns**

    


### fragmenstein.mpro.poised_pose_fx(pose)
Histidine in delta and cysteine in thiolate.


### fragmenstein.mpro.pose_fx(pose)
Histidine in delta.
