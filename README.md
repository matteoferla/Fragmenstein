# Fragmenstein
Scaffold hopping between bound compounds by stitching them together like a reanimated corpse.

## Premise

Aim: place a followup compound to the hits as faithfully as possible regardless of the screaming forcefields.
This makes a starting position for any subsequent calculations —say ∆∆G_bound vs. RMSD from the Frangmenstein pose after various optimisation cycles.

## Dramatis personae

There are **for now** two scripts here each with a namesake class.

* ``Fragmenstein`` makes the stitched together molecules
* ``Egor`` uses PyRosetta to minimise in the protein the fragmenstein followup.

## Description

Given a RDKit molecule and a series of hits it makes a spatially stitched together version of the initial molecule based on the hits.
Like Frankenstein's creation it violates the laws of chemistry. Planar trigonal topologies may be tetrahedral, bonds unnaturally long _etc._


* `.scaffold` is the combined version of the hits (rdkit.Chem.Mol object).
* `.chimera` is the combined version of the hits, but with differing atoms made to match the followup (rdkit.Chem.Mol object).
* `.positioned_mol` is the desired output (rdkit.Chem.Mol object)

Note, the hits have to be spatially aligned —i.e. extracted from crystal structures in bond form.

`.get_positional_mapping`, which works also as a class method, creates a dictionary of mol_A atom index to mol_B atom index
based on distance (cutoff 2&Aring;) and not MCS.

The code works in two broad steps, first a scaffold is made, which is the combination of the hits (by position).
Then the followup is placed. It is not embedded with constraint embed, which requires the reference molecule to have a valid geometry.
`.scaffold` and `.chimera` and `.positioned_mol` absolutely do not have this.
Novel side chains are added by aligning an optimised conformer against the closest 3-4 reference atoms.
Note that `.initial_mol` is not touched. `.positioned_mol` may have lost some custom properties, but the atom idices are the same.

## Algorithm

The following steps are done:

* `.merge_hits()`: merges the hits, the output `rdkit.Chem.Mol` object added as `.scaffold`.
* `.make_chimera()`: makes the atomic elements in `.scaffold` match those in the followup, the output `rdkit.Chem.Mol` object added as `.chimera`.
* `.place_followup()`  followup is places like the scaffold.

## Merger
The merger of the hits does not use MCS.
Instead all atoms are matched uniquely based on cartesian position with a cutoff of 2 Å (changeable).
The method ``get_positional_mapping`` works as a class method, so can be used for other stuff.
Here are the atoms in x0692 that map to x0305 and vice versa.

<img src="images/image0.svg" alt="x0692 common with x0305" width="200px">
<img src="images/image1.svg" alt="x0305 common with x0692" width="200px">

With maximum common substructure, a benzene ring can be mapped 6-ways, here none of that ambiguity is present although different issues arise.
Then the next step is fragment the second molecule to get the bits to add.

<img src="images/image2.svg" alt="fragments of x0305" width="200px">
<img src="images/image3.svg" alt="fragments of x0305" width="200px">
<img src="images/image6.svg" alt="fragments of x0305" width="200px">

The first fragment when added results in:

<img src="images/image3.svg" alt="fragments of x0305" width="200px">
<img src="images/image4.svg" alt="fragments of x0305" width="200px">
<img src="images/image5.svg" alt="fragments of x0305" width="200px">

Ditto for the second:

<img src="images/image8.svg" alt="fragments of x0305" width="200px">

And ditto again for a second fragment (x1249):

<img src="images/image14.svg" alt="fragments of x0305" width="200px">

## Elemental changes

During the elemental change, valence is taken into account resulting in appropriate positive charge.
This step is needed to avoid weird matches with the followup.

## Placing

To do a contrained embed in RDKit the reference need to have a good geometry. Consequently, this is not possible.
Therefore in the case of sidechains that are novel in the followup a optimised conformer is a aligned against the half placed followup
using the 3-4 atoms that are the closest neighbours within the half-placed structure and the side chain position copied from there for each bit.

<img src="images/grid.jpg" alt="fragments of x0305" width="400px">
<img src="images/overlay.png" alt="fragments of x0305" width="400px">

## Example

    hits = [Chem.MolFromMolFile(f'../Mpro/Mpro-{i}_0/Mpro-{i}_0.mol') for i in ('x0692', 'x0305', 'x1249')]
    followup = Chem.MolFromSmiles('CCNc1nc(CCS)c(C#N)cc1CN1C(CCS)CN(C(C)=O)CC1')
    monster = Fragmenstein(followup, hits)
    #monster = Fragmenstein(followup, hits, draw=True) for verbosity in a Jupyter notebook
    monster.make_pse('test.pse')
    
    display(monster.scaffold)
    display(monster.chimera) # merger of hits but with atoms made to match the to-be-aligned mol
    display(monster.positioned_mol) # followup aligned
    
    # further alignments... badly written way of doing this.
    monster.initial_mol = new_mol
    aligned = monster.place_followup(new_mol)
 
## Covalent

If the `Chem.Mol` has a dummy atom (element symbol: `*` within RDKit and smiles, but `R` in a mol file and PDB file) and
a molecule with a single atom is passed to `attachement` argument, then the covalent linker if absent in the hits is anchored
to that atom.
The default dummy atom can be overridden with `Fragmenstein.dummy:Chem.Mol` and `Fragmenstein.dummy_symbol:str`.

## Complicated MCS

Whereas the hit joining is done based on spatial overlaps. The followup is mapped to the blended scaffold by MCS.
First the list of possible MCS with really strict settings are found.
Then a set of mapping is sought which includes one of these by doing a new MCS search but very lax.
And going from very lax in increasing strictness. This prevents some weird mapping.

For more see `get_mcs_mapping`.
    
## Unresolved issues

Here is an example with a few issues.

<img src="images/unconnected.jpg" alt="unconnected" width="400px">

### Non-overlapping fragments

Non-overlapping fragments are discarded. Ideally they should be joined using the followup compound as a template _if possible_.
In the pictured case the SMILES submitted may not have been what was intended and should not be used connect the fragments.

### Imperfect projection

The projection approach is not perfect. The pictured example was affected by a bug (fixed), but this still is a problem in other cases.
This problem is quite apparent in the cases where atoms connecting to the sulfur are added:

![deviant](images/S_deviant.png)

The way the projection is done is via a single conformer.

### More than 4 templates

<img src="images/v_ugly.jpg" alt="ugly" width="400px">

When there are too many templates with large spread, these aren't merged resulting in a spiderweb scaffold.
This results in a non-unique mapping.


## Egor

_This script requires a module that I cannot share, but I was 70% through rewriting it, so should be rad._


Egor minimises the Fragmenstein monster in the protein using PyRosetta.

Egor has three minimisers that I tried:

* a modified cartesian FastRelax which works effectively.
* cartesian MinMover which gets stuck in local minimum in hard cases.
* PertMinMover which behaves weirdly...

<img src="images/movers.jpg" alt="movers" width="400px">

* The template ought to be relaxed beforehand —but this can be skipped.
* An ideal conformer set made into a `params` file.

Both of which ATM are done with a different script.


    e = Egor(pose, constraint_filename)
    e.minimise(10)

Where pose is a `pyrosetta.Pose` instance.
But Egor can be initialised with `Egor.from_pdbfile(..)` or `Egor.from_pdbblock(..)`.
The latter is nothing more than:

    e.repack_neighbors()
    mover = e.get_FastRelax(10)
    # mover = e.get_PertMinMover()
    # mover = e.get_MinMover()
    mover.apply(e.pose)

## See also

* [my messy code for Covid19 moonshot](https://github.com/matteoferla/SARS-CoV-2_CL3_covalent_docking).

