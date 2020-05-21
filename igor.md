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
The latter is nothing more than:

    e.repack_neighbors()
    mover = e.get_FastRelax(10)
    # mover = e.get_PertMinMover()
    # mover = e.get_MinMover()
    mover.apply(e.pose)

Do note that when Victor calls Igor the constraint file used will have the atoms that should be constrained by position.
The `Igor.coordinate_constraint` (default = 1) controls the stringency.
Note that the `get_mod_FastRelax` has a weight in case it runs of constrain to position alone.
The coordinate constraint is a harmonic function of standard deviation 1 plus the atom `_StDev` property
—total math fudge: I am happy to hear better suggestions!