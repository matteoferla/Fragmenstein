## Results reproducibility.

In order to obtain reproducible monster results you can set the monster random seed:

```
Monster(hits=hits, random_seed=132)
```

In the same manner, to obtain variability, just try different seeds.

For Victor, in order to ensure reproducibility, the random seed needs to be fixed at the
`pyrosetta.init` call and the `monster_random_seed` argument needs to be also set.

```
pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false -constant_seed true -jran 987') 
    Victor(hits=hits, pdb_filename=pdb_filename, monster_random_seed=987)
```