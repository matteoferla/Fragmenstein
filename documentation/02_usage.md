# Pythonic usage

### Classes

There are four main classes —named after characters from the Fragmenstein book and movies:

* `Monster` makes the stitched together molecules indepently of the protein — [documentation](monster/monster.md)
* `Igor` uses PyRosetta to minimise in the protein the fragmenstein monster followup — [documentation](further-detail/igor.md)
* `Victor` is a pipeline that calls the parts, with several features, such as warhead switching —[documentation](further-detail/victor.md)
* `Laboratory` does all the combinatorial operations with Victor (specific case)

NB. In the absence of `pyrosetta` (which requires an academic licence), all bar ``Igor`` work and 
alternative Victor classes need to be used, for example 
`Wictor` (RDkit minimisation only), `OpenVictor (using OpenMM).

Additionally, there are a few minor classes.

One of these is ``mRMSD``, a multiple RMSD variant which does not superpose/align and bases which atoms 
to use on coordinates —[documentation](further-detail/mrmsd.md)

The class `Walton` performs geometric manipulations of compounds, to set them up to demonstrate
features of Fragmenstein (like captain Walton, it does not partake in the plot, but is key to the narration)

There are two module hosted elsewhere:

* ``Rectifier`` from [molecular_rectifier](https://github.com/matteoferla/molecular_rectifier) is a class that corrects mistakes in the molecule automatically merged by ``Monster``.
* ``Params`` from [rdkit to params module](https://github.com/matteoferla/rdkit_to_params) parameterises the ligands

### Combine
It can also merge and link fragment hits by itself and find the best scoring mergers.
For details about linking see [linking notes](further-detail/linking.md).
It uses the same overlapping position clustering, but also has a decent amount of impossible/uncommon chemistry prevention.

Monster:

```python
from fragmenstein import Monster
hits: List[Chem.Mol] = ... # 1 or more RDKit.Chem.Mol, sanitised, w/ conformer, preferably without explicit Hs
monster = Monster(hits=hits)
monster.combine()
monster.positioned_mol #: RDKit.Chem.Mol
```

Victor:

```python
from fragmenstein import Victor
import pyrosetta
pyrosetta.init( extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

hits: List[Chem.Mol] = ...
victor = Victor(hits=hits, # List of 1 or more RDKit.Chem.Mol
            pdb_filename='foo.pdb',  # or pdb_block='ATOM 1 MET ...'
            covalent_resi=1) # if not covalent, just put the first residue or something.
victor.combine()
victor.minimized_mol
```
The PyRosetta init step can be done with the helper function:
```python
Igor.init_pyrosetta()
```

The two seem similar, but Victor places with Monster and minimises with Igor.
As a result it has ∆G_bind energy score (difference between holo minus apo+ligand Gibbs free energy predictions):

victor.ddG

Fragmenstein is not really a docking algorithm as it does not find the pose with the **lowest energy** 
within a given volume.
Consequently, it is a method to find how **faithful** is a given followup to the hits provided.
Hence the minimised pose should be assessed by the RMSD metric or similar
and the ∆G_bind score used solely as a cutoff —lower than zero.

For a large number of combination:

```python
from fragmenstein import Laboratory

pdbblock: str = ... # a PDB block
hits: List[Chem.Mol] = ... # 1 or more RDKit.Chem.Mol
lab = Laboratory(pdbblock=pdbblock, covalent_resi=None)
combinations:pd.DataFrame = lab.combine(hits, n_cores=28)
```

### Place
Here is [an interactive example of placed molecules](https://michelanglo.sgc.ox.ac.uk/r/fragmenstein).

It is rather tolerant to erroneous/excessive submissions (by automatically excluding them)
and can energy minimise strained conformations.
![summary](../images/new_summary.jpg)

Three mapping approaches were tested, but the key is that hits are pairwise mapped to each other by means 
of one-to-one atom matching based upon position as opposed to similarity which is easily led astray. 
For example, note here that the benzene and the pyridine rings overlap, not the two pyridine rings:

<img src="../images/position_over_mcs.jpg" width="300px">

### RDkit only and OpenMM

PyRosetta is needed for the pocket-centric minimisation.
Two alternatives are available:

* `Wictor` (without): stops at the RDKit minimisation
* `OpenVictor` (with OpenMM): uses OpenMM to minimise in the protein

Whereas the PyRosetta steps operate via Igor, OpenVictor uses Fritz.
OpenMM is a lot slower than PyRosetta on CPU only,
but is free, open source and potentially more accurate.

Igor is a much larger class as it needs to disable rotamer sampling and other things,
which is not an issue in OpenMM.

A further detail is that openMM is already parallel,
therefore when using with `Laboratory` request only one core.
```python
from fragmenstein import Laboratory, OpenVictor
Laboratory.Victor = OpenVictor

lab = Laboratory(pdbblock=MPro.get_template())
combinations: pd.DataFrame = lab.combine(hits,
                                     n_cores=1,  # 1 core unless $OPENMM_CPU_THREADS is set
                                     timeout=600,  # 2 minutes
                                     combination_size=2,  # pairwise
                                     max_tasks=0)  # 0 is no chunking
```

## Subclassing

`Laboratory` uses `Victor`, but uses a class attribute `.Victor` pointing to the `Victor` class,
when running. Likewise, `Victor` has `.Monster`. There is no `Igor` equivalent as the Igor calls by Victor 
are far from generic.

There are some empty methods aimed at easier subclassing:

* `Monster.post_ff_addition_step` —called after the MMFF forcefield is added in the `Monster.mmff_minimize` method
* `Victor.post_monster_step` -called in the `Victor._calculate_*_chem` methods after the `Monster` is stitched together
* `Victor.post_params_step` - called in the `Victor._calculate_*_thermo` methods after the `Params` are added
* `Victor.pre_igor_step` - called in the `Victor._calculate_*_thermo` methods before the `Igor` setup is started
* `Victor.pose_mod_step` - called in the `Victor._calculate_*_thermo` methods after the pose in loaded, an alternative to `pose_fx`
* `Victor.post_igor_step` - called in the `Victor._calculate_*_thermo` methods after the `Igor` minimisation is finished

## E-amide isomerism

The static method `Monster.inspect_amide_torsions` called on a Chem.Mol will return a list of the amide torsions.
This is because some mergers do result in cis (E) amides, which are not ideal.
By default, the MMFF forcefield will penalise exocyclic cis-amides.
There is a setting, `ff_prevent_cis` (as env `$FRAGMENSTEIN_FF_PREVENT_CIS`),
which when set to `false` will disable the addition of a penalty against cis-amides to the MMFF forcefield
when `Monster.mmff_minimize` is called.