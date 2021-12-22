# Victor

Victor is the overarching class. This has many features and can be rather complicated in its setup: but then a point and click
solution that works universally without customisation is a bad solution for molecular modelling.

## Combine vs. Place

Like Monster Victor has two main calculation methods: `combine` and `place`.
Whereas in Monster these are independent of protein neighbourhood,
in Victor they are not thanks to Igor. 

Victor does the following steps:

* Is given a followup to test and the mol objects of its inspiration and the pdb template file.
* Calls Monster class
* Parameterises the mol
* Generates the constraints
* Calls Igor
* Can do some extras

## Class attributes
Core settings controlling its behaviour can be set via class attribute. For example,
the warhead conversion methods can be controlled via:

* `covalent_definitions` (currently defined only for cysteine)
* `warhead_definitions`

While other include:

* `error_to_catch`
* `constraint_function_type`
* `work_path`
    
## Class methods

Two key class methods are `Victor.extract_mols` and `Victor.extract_mol`.
Also several methods for covalent operations — see [covalent notes](covalents.md)

Note that Igor has pyrosetta specific class methods, e.g. downloading electron density for template prep etc.

## Example

> The following code may have changed. And several of these methods are

Here is a real world usage that uses multiple features:

Import pyrosetta and initialised before everything (applies to Igor too):

    import pyrosetta
    pyrosetta.init(extra_options='-no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false')
    from fragmenstein import Igor, Monster, Victor
    import logging, csv, json
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
configure whither to save and whence to load:

    Victor.work_path = '../Mpro_fragmenstein'
    Victor.enable_stdout(logging.WARNING)
    mpro_folder = '/Users/matteo/Coding/rosettaOps/Mpro'  
        
alternatively `Victor.enable_logfile('reanimate.log',logging.DEBUG)` (see [logging notes](logging_and_debugging.md)).
    
add extra constraints that are warhead & protein specific.
note that the warhead definitions contain preferred names for the connecting atoms and their neighbours

    for cname, con in [('chloroacetamide', 'AtomPair  H  145A  OY  1B HARMONIC 2.1 0.2\n'),
                      ('nitrile', 'AtomPair  H  145A  NX  1B HARMONIC 2.1 0.2\n'),
                      ('acrylamide', 'AtomPair  H  143A  OZ  1B HARMONIC 2.1 0.2\n'),
                      ('vinylsulfonamide', 'AtomPair  H  143A  OZ1 1B HARMONIC 2.1 0.2\n')
                      ]:
        Victor.add_constraint_to_warhead(name=cname, constraint=con)
        
Here is the definition of a nitrile warhead, for example:

    {'name': 'nitrile',
    'covalent': 'C(=N)*',  # zeroth atom is attached to the rest
    'covalent_atomnames': ['CX', 'NX', 'CONN1'],
    'noncovalent': 'C(#N)',  # zeroth atom is attached to the rest
    'noncovalent_atomnames': ['CX', 'NX']
    }
    
This allows warheads to be mixed and matched.
<img src="images/warheads.jpg" alt="warhead" width="400px">

The choice of the protein template is a bit weak.
I plan to experiment with minimisation against averaged electron densities.
    
    def get_best(hit_codes):
        return Victor.closest_hit(pdb_filenames=[f'{mpro_folder}/Mpro-{i}_0/Mpro-{i}_0_bound.pdb' for i in hit_codes],
                            target_resi=145,
                            target_chain='A',
                            target_atomname='SG',
                            ligand_resn='LIG')
    
There is a change I require to the pose
    
    def pose_fx(pose):
            pose2pdb = pose.pdb_info().pdb2pose
            r = pose2pdb(res=41, chain='A')
            MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
            MutateResidue(target=r, new_res='HIS').apply(pose)
    
Define all the steps
    
    def reanimate(smiles, name, hit_codes):
        hits = [get_mol(i) for i in hit_codes]
        best_hit = get_best(hit_codes)
        Victor.journal.debug(f'{name} - best hit as starting is {best_hit}')
        apo = best_hit.replace('_bound', '_apo-desolv')
        print(f'reanimate(smiles="{smiles}", name="{name}", hit_codes={hit_codes})')
        reanimator = Victor(hits=hits,
                            pdb_filename=apo,
                            ligand_resn='LIG',
                            ligand_resi='1B',
                            covalent_resn='CYS', covalent_resi='145A',
                            pose_fx = pose_fx
                            )
        reanimator.place(smiles=smiles,
                         extra_constraint='AtomPair  SG  145A  NE2  41A HARMONIC 3.5 0.2\n',
                         long_name=name)
        return reanimator
     
Read the data and do all warhead combinations if covalent. This data is actually from
[github.com/postera-ai/COVID_moonshot_submissions](https://github.com/postera-ai/COVID_moonshot_submissions).
For which there is a method in `fragmenstein.mpro`.
       
    data = csv.DictReader(open('../COVID_moonshot_submissions/covid_submissions_all_info.csv'))
    
    issue = []
    
    for row in data:
        if row['covalent_warhead'] == 'False':
            pass
            reanimate(name = row['CID'], hit_codes = row['fragments'].split(','), smiles=row['SMILES'])
        else:
            print(f'Covalent: {row["CID"]}')
            for category in ('acrylamide', 'chloroacetamide', 'vinylsulfonamide', 'nitrile'):
                if row[category] == 'True':
                    combinations = Victor.make_all_warhead_combinations(row['SMILES'], category)
                    if combinations is None:
                        issue.append(row["CID"])
                        break
                    for c in combinations:
                        reanimate(name = row['CID']+'-'+c, hit_codes = row['fragments'].split(','), smiles=combinations[c])
                    break
            else:
                print(f'What is {row["CID"]}')
                issue.append(row["CID"])
      
## Methods for inheritance
The above could have been customised further, by making a class that inherits Victor and defining
 `post_params_step`, `post_fragmenstein_step`, `pose_mod_step` or `post_igor_step`, which are empty methods
intended to make subclassing Victor easier as these are meant to be overridden
—NB `pose_mod_step` is run if not `pose_fx` is given.

## Coordinate constraint
The coordinate constraint generated for Igor, the minimiser can be changed from `HARMONIC` (x = mean)
to `FLAT_HARMONIC` (tol = max distance of contributing atoms) and `BOUNDED` (fixed penalty potential well).

Do note that the reference atom for the constraint is the covalent residue, regardless of whether the ligand is covalent.
