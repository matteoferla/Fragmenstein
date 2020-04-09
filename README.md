# Fragmenstein
Scaffold hopping between bound compounds by stitching them together like a reanimated corpse.

## Premise

Aim: place a followup compound to the hits as faithfully as possible regardless of the screaming forcefields.

## Description

Given a RDKit molecule and a series of hits it makes a spatially stitched together version of the initial molecule based on the hits.


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

## Example

    hits = [Chem.MolFromMolFile(f'../Mpro/Mpro-{i}_0/Mpro-{i}_0.mol') for i in ('x0692', 'x0305', 'x1249')]
    probe = Chem.MolFromSmiles('CCNc1nc(CCS)c(C#N)cc1CN1C(CCS)CN(C(C)=O)CC1')
    monster = Fragmenstein(probe, hits)
    #monster = Fragmenstein(probe, hits, draw=True) for verbosity in a Jupyter notebook
    monster.make_pse('test.pse')
    
    display(monster.scaffold)
    display(monster.chimera) # merger of hits but with atoms made to match the to-be-aligned mol
    display(monster.positioned_mol) # followup aligned
    
    # further alignments... badly written way of doing this.
    monster.initial_mol = new_mol
    aligned = monster.place_probe(new_mol)

For details

## See also

* [my messy code for Covid19 moonshot](https://github.com/matteoferla/SARS-CoV-2_CL3_covalent_docking).

