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
                           regex_name='x(\d+)', # optional regex to run get rid of fluff
                           proximityBonding=False, # change this to True if you lack CONECT entries (disconnected ligand)
                           )
```

To have a gander
```jupyterpython
from rdkit import Chem
Chem.Draw.MolsToGridImage(mols.values())
```
                           
NB. there is a Victor method called `from_files`, 
this is for restarting a run precariously from the saves and has nothing to do with these.

NB2. One thing to be vigilant for is that the ligands are in aligned protein.
Fragmenstein does nothing to ensure this is true as multichain protein etc. make everything harder and
it does not require much effort at all to fix (see Troubleshooting).

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

### Troubleshooting
To see the Pyrosetta pose
```jupyterpython
import nglview as nv

view = nv.show_rosetta(victor.igor.pose)
view
```
If there is something like this:
![farway](../images/farway.png)
Then it means that the PDB where the molecule was extracted is in a different from the template.
To fix, before extracting or before using a template align them in say PyMOL (`align` command).
Just remember than PyMOL strips LINK entries making covalents non-covalent
(unless `proximityBonding=True` is used in the extraction, which is not recommended).
Doing it within python:
```jupyterpython
import pymol2

with pymol2.PyMOL() as pymol:
    pymol.cmd.load('minimized.pdb', 'mini')
    pymol.cmd.load('x01234_refined.pdb', 'hit')
    pymol.cmd.align('mini', 'hit and chain A')
    pymol.cmd.save('moved.pdb', 'mini')
```

Another issue may arise when hydrogens are present in the hits somehow.

If the calculations fail along a step the compounds can be inspected, but these are not `victor.minimised_mol`

If an error happens polishing up the minimised molecule from pyrosetta:

```jupyterpython
ligand = victor.igor.mol_from_pose()
from rdkit import Chem
from rdkit.Chem import AllChem
AllChem.Compute2DCoords(ligand)
ligand = AllChem.RemoveAllHs(ligand)
ligand
```

If the error happens during Igor, but Monster worked fine, the ligand is in `victor.monster.positioned_mol`.
Additionally if the issue is with one of the atoms and the indices are required there is a debug focused method in Monster:

```jupyterpython
mol = victor.monster.positioned_mol
victor.monster.draw_nicely(mol)
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

## Michelanglo
To make an interactive page in [Michelanglo](https://michelanglo.sgc.ox.ac.uk/), like [this example](https://michelanglo.sgc.ox.ac.uk/r/fragmenstein).
one can use the [michelanglo_api](https://github.com/matteoferla/MichelaNGLo-api) (pip name is `michelanglo-api`).
The data is stored in a github repo. For a detailed example see [pipeline](pipeline.md).

First, make an SDF of the results with properties to show in the table say via `rdkit.Chem.PandasTools`.
This requires a first compound with metadata values. This will change one day.

then make or get a Michelanglo page
```jupyterpython
from michelanglo_api import MikeAPI
mike = MikeAPI('username', 'password')
page = mike.convert_pdb('6WOJ') # make new
# or...
#p = mike.get_page('xxxxxxxx')  # retrieve old
page.retrieve()
page.show_link()
```
Fix up... etc.
```jupyterpython
page.description = 'Hello world. '
page.loadfun = ''
page.columns_viewport = 6
page.columns_text = 6
```
Add the data
```jupyterpython

gitfolder='/Users/you/path_to_your_github_repo_on_your_machine'
sdfile='/Users/you/path_to_sdfile.sdf'
folder = 'folder_name_within_repo'
targetfolder=f'{gitfolder}/{folder}'

# move the sdf_file to individual mol files in your repo
page.sdf_to_mols(sdfile=sdfile,
             targetfolder=targetfolder,
             skip_first=True) # first row is metadata in a SDF for XChem
# make a json for the table
page.sdf_to_json(sdfile=sdfile,
             keys=('∆∆G', 'comRMSD', 'N_constrained_atoms', 'runtime', 'disregarded', 'smiles'),
             key_defaults=(999., 999., 0, 999., 'NA', 'NA'), #what to set stuff that is null
             filename=f'{targetfolder}/data.json')
# make a table
page.make_fragment_table(sdfile=sdfile,
               username='matteoferla',
               repo_name='Data_for_own_Michelanglo_pages',
               foldername=folder,
               protein_sele='145:A', # show this on protein. NGL selection
               sort_col=2, #sort by column index 2.
               sort_dir='asc', #asc or desc
               template_row=-1, # is the template a file called `template.pdb` (-1) or a filename in the row n?
               fragment_row=1, # the inspiration fragments (-1 for none). The names must match with or without a .mol.
               jsonfile='data.json')
# commit changes
page.commit()
```



