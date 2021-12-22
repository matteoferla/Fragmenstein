# Fragmenstein
Scaffold hopping between bound compounds by stitching them together like a reanimated corpse.
<img src="images/fragmenstein.jpg" width="300px">


## Stitched molecules

Fragmenstein can perform two different tasks.

* **Combine** hits
* **Place** a given followup molecule (SMILES) based on series of hits

Like Frankenstein's creation it may violate the laws of chemistry.
Planar trigonal topologies may be tetrahedral, bonds unnaturally long _etc._
This monstrosity is therefore then energy minimised with strong constraints within the protein.

## Classes

There are three main classes, named after characters from the Fragmenstein book and movies:

* ``Monster`` makes the stitched together molecules indepent of the protein — [documentation](documentation/monster/monster.md)
* ``Igor`` uses PyRosetta to minimise in the protein the fragmenstein followup — [documentation](documentation/igor.md)
* ``Victor`` is a pipeline that calls the parts, with several features, such as warhead switching —[documentation](documentation/victor.md)

NB. In the absence of `pyrosetta` (which requires an academic licence), all bar ``Igor`` work.

Additionally, there are a few minor classes.

One of these is ``mRMSD``, a multiple RMSD variant which does not align and bases which atoms 
to use on coordinates —[documentation](documentation/mrmsd.md)

There are two module hosted elsewhere:

* ``Rectifier`` from [molecular_rectifier](https://github.com/matteoferla/molecular_rectifier) is a class that corrects mistakes in the molecule automatically merged by ``Monster``.
* ``Params`` from [rdkit to params module](https://github.com/matteoferla/rdkit_to_params) parameterises the ligands

### Combine
It can also merge and link fragment hits by itself and find the best scoring mergers.
For details about linking see [linking notes](documentation/linking.md).
It uses the same overlapping position clustering, but also has a decent amount of impossible/uncommon chemistry prevention.

Monster:
 
    from fragmenstein import Monster
    monster = Monster(hits=[hits_a, hit_b])
    monster.combine()
    monster.positioned_mol # RDKit.Chem.Mol
    
Victor:

    from fragmenstein import Victor
    import pyrosetta
    pyrosetta.init( extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

    victor = Victor(hits=[hits_a, hit_b], 
                    pdb_filename='foo.pdb',  # or pdb_block='ATOM 1 MET ...'
                    covalent_resi=1) # if not covalent, just put the first residue or something.
    victor.combine()
    victor.minimized_mol
    
The two seem similar, but Victor places with Monster and minimises with Igor.
As a result it has energy scores

    victor.ddG
    
Fragmenstein is not really a docking algorithm as it does not find the pose with the **lowest energy** 
within a given volume.
Consequently, it is a method to find how **faithful** is a given followup to the hits provided.
Hence the minimised pose should be assessed by the RMSD metric or similar
and the ∆∆G score used solely as a cutoff —lower than zero.    

## Place
Here is [an interactive example of placed molecules](https://michelanglo.sgc.ox.ac.uk/r/fragmenstein).

It is rather tolerant to erroneous/excessive submissions (by automatically excluding them)
and can energy minimise strained conformations.
![summary](images/new_summary.jpg)

Three mapping approaches were tested, but the key is that hits are pairwise mapped to each other by means 
of one-to-one atom matching based upon position as opposed to similarity which is easily led astray. 
For example, note here that the benzene and the pyridine rings overlap, not the two pyridine rings:

<img src="images/position_over_mcs.jpg" width="300px">

### Examples

Monster:
 
    from fragmenstein import Monster
    monster = Monster(hits=[hits_a, hit_b])
    monster.place_smiles('CCO')
    monster.positioned_mol
    
Victor:

    from fragmenstein import Victor
    import pyrosetta
    pyrosetta.init( extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

    victor = Victor(hits=[hits_a, hit_b], pdb_filename='foo.pdb')
    victor.place('CCO')
    victor.minimised_mol
    
For a lengthier example see [example notes](documentation/example.md).

### MPro example
To use SAR-COV-2 MPro as a test bed, the following may be helpful:

* `fragmenstein.MProVictor`, a derived class (of `Victor`), with various presents specific for MPro.
* `fragemenstein.get_mpro_template()`, returns the PDB block (str) of MPro
* `fragemenstein.get_mpro_molblock(xnumber)`, returns the mol block (str) of a MPro hit from Fragalysis
* `fragemenstein.get_mpro_mol(xnumber)`, as above but returns a `Chem.Mol` instance.

## Other features

* [Covalent hits](documentation/covalents.md)
* [Logging](documentation/logging_and_debugging.md)

## Installation

### Requires RDKit
To install for system Python3 on Linux:

    sudo apt-get install python3-rdkit librdkit1 rdkit-data
    
To install for system Python3 on MacOS via Brew:

    brew install rdkit --with-python3
    
To install for Conda Python3

    conda install -c conda-forge rdkit
    
### Requires Pyrosetta

Pyrosetta requires a password to be downloaded (acamedic licence) obtained by https://els2.comotion.uw.edu/product/pyrosetta. This is a different licence from the Rosetta one. The username of the latter is formatted variant of "academic user", while the former is the name of a researcher whose names bares an important concept in protein folding. Pyrosetta can be downloaded via a browser from http://www.pyrosetta.org/dow. Or with a variant of the following:

    
    curl -u 👾👾👾:👾👾👾https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python38.linux/PyRosetta4.Release.python38.linux.release-273.tar.bz2 -o a.tar.bz2
    tar -xf a.tar.bz2
    cd PyRosetta4.Release.python38.linux
    sudo pip3 install .

### Fragmenstein and dependencies

Install from pipy

    sudo pip3 install fragmenstein

## Origin

> See [Fragmenstein and COVID moonshot](documentation/covid.md).

Fragmenstein was created to see how reasonable are the molecules of fragment mergers submitted
in [the COVID moonshot project](https://discuss.postera.ai/c/covid), because after all the underlying method is 
fragment based screening.
[This dataset](https://github.com/postera-ai/COVID_moonshot_submissions) has some unique peculiarities that potentially
are not encountered in other projects.

## Autogenerated documentations

For more see the source code or the [Sphinx converted documentation](documentation/sphinx-docs.md).

## Changes

Some changes to the algorithm may happen, 
see [changelog](documentation/changelog_0.6.md) and [wip.md](documentation/wip.md) for more.


