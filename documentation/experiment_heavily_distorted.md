## Experiment: heavily distorted molecule

To test the system I made a really nastily distorted xylene and 
checked what would happen passing it to Fragmenstein.

Make a distorted xylene:

```python
from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.Chem import AllChem
import random
from rdkit_to_params import Params

xylene = Chem.MolFromSmiles('c1ccc(C)cc1C')
AllChem.EmbedMolecule(xylene)
Params.load_mol(xylene, name='XYL')  # adds atomnames in place basically
Chem.MolToMolFile(xylene, 'xylene.pre.mol')
conf = xylene.GetConformers()[0]
alpha, beta = 0.3, 0.5 # scale & shape
for ai in range(xylene.GetNumAtoms()): #type: int
    point :Point3D = conf.GetAtomPosition(ai)
    point.x += random.weibullvariate(alpha,beta)
    point.y += random.weibullvariate(alpha,beta)
    point.z += random.weibullvariate(alpha,beta)
    conf.SetAtomPosition(ai, point)
Chem.MolToMolFile(xylene, 'xylene.post.mol')
```

Now, let's make a PDB file with a calcium 100 Å from the origin.
This is because the coordinate constraints require a second point
in order to fix to the origin basically.
As a result, the Fragmenstein protocol cannot work with out a template.

```python
import pyrosetta
import pyrosetta_help as ph

pyrosetta.init(extra_options=ph.make_option_string(no_optH=False,
                                                   ignore_unrecognized_res=True,
                                                   load_PDB_components=False,
                                                   ignore_waters=False),
               set_logging_handler='logging'
               )
pose = pyrosetta.pose_from_sequence('X[CA]')
res = pose.residue(1)
far = pyrosetta.rosetta.numeric.xyzVector_double_t(100,0,0)
for i in range(1, 1+res.natoms()):
    xyz = res.xyz(i)
    xyz.x += 100
    res.set_xyz(i, xyz)
pose.dump_pdb('ref.pdb')
```

Now lets feed it the janky xylene

```python
from fragmenstein import Victor
import nglview as nv

victor = Victor(hits=[xylene,],
                pdb_filename='ref.pdb',
                ligand_resn='XYL',
                covalent_resi=1)
# victor.monster_mmff_minisation = False/True
victor.combine()
# show it!
view = nv.show_file(StringIO(victor.unminimised_pdbblock), ext='pdb')
view.update_ball_and_stick(colorValue='#00B4C4', smoothSheet=True)  # mid-dark cyan
c2 = view.add_component(StringIO(victor.minimised_pdbblock), ext='pdb')
c2.update_ball_and_stick(colorValue='#F8766D', smoothSheet=True)   # salmon
view
```
Whereas `place` stretches into place a provided idealised molecule,
`combine` stretches out of place into a more ideal state a fragmenstein molecule.
So the combined janky xylene with nothing else will be ugly.
To speed things up and to work without PyRosetta
an optional pre-minimisation step using rdkit is done which
can be disabled via `victor.monster_mmff_minisation = False`.

If PyRosetta is run, the result is the same as the FastRelax call
uses `ref2015_cart` scorefunction with cartersian extra enabled.

If the initial distorsion is too much it fails as 4 Å is too long.