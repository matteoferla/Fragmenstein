# Examples

## Config
Victor is the main entrypoint to the module.

```jupyterpython
from fragmenstein import Victor
```

### Logging
set logging to stdout
```jupyterpython
import logging
Victor.enable_stdout(logging.INFO)
```

or to file
```jupyterpython
import logging
Victor.enable_logfile('test.log', logging.INFO)
```

the Logger instance is `Victor.journal`, and it captures rdkit and pyrosetta logs 
if `enable_stdout` or `enable_logfile` are called with the argument `capture` set to `True` (default)
or the method `capture_logs` is called directly.

### Pyrosetta

Pyrosetta needs to be initialised as normal

```jupyterpython
import pyrosetta
pyrosetta.init(extra_options='-no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false')
```

Alternatively, a cleaner way can be using a helper function in laboratory

```jupyterpython
import pyrosetta
from fragmenstein.laboratory import make_option_string
extras = make_option_string(no_optH=False,
                            mute='all',
                            ignore_unrecognized_res=True,
                            load_PDB_components=False)
pyrosetta.init(extra_options=extras)
```

If you want to just use `Monster` and do not have `pyrosetta` installed for licencing reasons,
the following works fine.

```jupyterpython
from fragmenstein import Monster
```

If you have pyrosetta installed, but want to mimic this behaviour (??)

```jupyterpython
import sys
sys.modules['pyrosetta'] = None
# this raises on `import pyrosetta` a ModuleNotFoundError
```
### Output path

set working path
```jupyterpython
Victor.work_path = 'test'
```
### Hits

Create mol files of reference hits
given a single PDB file return the RDKit Chem.Mol molecule with covalent as `*`

```jupyterpython
mol = Victor.extract_mol(name='x01234',
                         filepath='here.pdb',
                         smiles='CCO',
                         ligand_resn='LIG',
                         removeHs=False)
```

To parse a whole folder (flat) of PDBs and get a dictionary of names -> Chem.Mol
```jupyterpython
mols = Victor.extract_mols(folder='PDBs',
                           smilesdex={'x01234': 'CCO'}, # optional
                           ligand_resn= 'LIG',
                           regex_name='x(\d+)', # optional regex to run get rid of fluff)
```
                           
NB. there is a Victor method called `from_files`, 
this is for restarting a run precariously from the saves and has nothing to do with these.

## Template
Like in a docking experiment, the template protein conformer is important
even if sidechains and backbones can move.
To minimise with Pyrosetta, it is best to minimise a ligand bond protein and remove the ligand
afterwards.

The ligand needs to be parameterised first. Example
(for more see [rdkit_to_params](https://github.com/matteoferla/rdkit_to_params))

```jupyterpython
from rdkit_to_params import Params

params = Params.from_smiles_w_pdbfile(pdb_file='5BV6_clean.pdb',
                              smiles='c1nc2c(n1[C@H]3[C@@H]([C@H]4[C@H](O3)CO[P@@](=O)(O4)[O-])O)N=C(NC2=O)N',
                            name='35G')
params.dump('35G.params')
```
Have a gander to see all is good
```jupyterpython
import nglview
nglview.show_rosetta(params.test())
```
Load pose:
```jupyterpython
from typing import *


def get_pose(pdb_filename: str,
             params_filenames: Optional[List[str]]) -> pyrosetta.Pose:
    """A fxn to load a pose, with a list of params"""
    pose = pyrosetta.Pose()
    if params_filenames and isinstance(params_filenames, pyrosetta.rosetta.utility.vector1_string):
        pyrosetta.generate_nonstandard_residue_set(pose, params_filenames)
    if params_filenames and isinstance(params_filenames, list):
        params_filenames2 = pyrosetta.rosetta.utility.vector1_string()
        params_filenames2.extend(params_filenames)
        pyrosetta.generate_nonstandard_residue_set(pose, params_filenames2)
    else:
        pass
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdb_filename)
    return pose


pose = get_pose(pdb_filename='5BV6_clean.pdb',
             params_filenames=['35G.params'])
```
Igor has a function to relax against ED (taken from [here](http://blog.matteoferla.com/2020/04/how-to-set-up-electron-density.html))

```jupyterpython
from fragmenstein import Igor
Igor.relax_with_ED(pose, 'map.ccp4')
pose.dump_pdb('relaxed.pdb')
pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose)
pose.dump_pdb('template.pdb')
```
## Merge

Victor has two modes, merge and place.
Say we extracted `hit_a` and `hit_b` as `Chem.Mol` instances:

```jupyterpython
victor = Victor(hits=[hits_a, hit_b], pdb_filename='template.pdb')
# victor.place('CCO') # to place.
victor.combine()
victor.minimised_mol
```
Using code from [here](http://blog.matteoferla.com/2021/02/multiple-poses-in-nglview.html)
one could inspect the starting hits in NGLView in the notebook etc.

Alternatively,
```jupyterpython
victor.make_pse()
```
One could see the numbers with:
```jupyterpython
victor.summarize()
```

Or individual values

```jupyterpython
victor.ddG
```

## Advanced

### Warheads

See [covalents](covalents.md)

### Monster

A problem is that the combined molecules may not be enamine purchasable.
So the placed molecule could be a purchaseable.
So getting a positioned mol

```jupyterpython
from fragmenstein import Monster
from rdkit import Chem

monster = Monster(hits=[hits_a, hit_b])
monster.combine()
smiles = Chem.MolToSmiles(monster.positioned_mol)
```

Using [enamine-real-search API](https://github.com/xchem/enamine-real-search)

```jupyterpython
from search import EnamineSession

session = EnamineSession()
similarity_results = session.similarity_search(smiles=smiles, threshold=0.1)
```

```jupyterpython

victor = Victor(hits=[hits_a, hit_b], pdb_filename='template.pdb')
victor.place(similarity_results.smiles) # to place.
```




