# Igor

Igor minimises the Fragmenstein monster in the protein using PyRosetta.

Igor uses the package [rdkit_to_params](https://github.com/matteoferla/rdkit_to_params) to parameterise the compounds.

Igor has three minimisers that I tried:

* a modified cartesian FastRelax which works effectively.
* cartesian MinMover which gets stuck in local minimum in hard cases.
* PertMinMover which behaves weirdly...

<img src="images/movers.jpg" alt="movers" width="400px">

Igor gets run as follows:

    e = Igor(pose, constraint_filename)
    e.minimise(10)

Where pose is a `pyrosetta.Pose` instance.
But Igor can be initialised with `Igor.from_pdbfile(..)` or `Igor.from_pdbblock(..)`.

Do note that when Victor calls Igor the constraint file used will have the atoms that should be constrained by position.
The `Igor.coordinate_constraint` (default = 1) controls the stringency.
Victor will make Igor repeat the minimisation with increasingly lower stringency until the ∆∆G is negative 
or when the weight is near zero.
In the absence of coordinate constraints in the constraint file (Victor generates it from the constraints),
the argument `default_coord_constraint=True` passed to `.minimise` will set `relax.constrain_relax_to_start_coords(True)`.
Note that the `get_mod_FastRelax` has a weight in case it runs of constrain to position alone.
As Victor sets the constraints, the class attribute of Victor `constraint_function_type` can be `HARMONIC`, `FLAT_HARMONIC` (default) or `FADE`.
