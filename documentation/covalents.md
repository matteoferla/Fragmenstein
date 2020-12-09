## Covalent

Covalent forms of a molecule are marked within the `Chem.Mol` instances in fragmenstein with a "dummy atom".
This is element symbol `*` within RDKit and smiles, but `R` in a mol file and PDB file.

This is essential for the params file (topology file for the ligand).

Consequently, hits are best extracted with `Victor.extract_mols` or `Victor.extract_mol`.
Fragalysis does not do this, so the hits from there will be treated as regular compounds.

## Monster

A molecule with a single atom is passed to `attachement` argument of `Monster`,
then the covalent linker if absent in the hits is anchored to that atom.
If you know a ligand that is covalent `Victor.find_attachment` can be used to extract the attachment atom.

## Warheads

Victor has a cysteine reactive warheads operations
In the class attribute `.warhead_definitions` are stored conversions and atom names
(cf. MProVictor for an example).

**acrylamide** `C(=O)C=C` => `C(=O)CC*`
**chloroacetamide** `C(=O)C[Cl]` => `C(=O)C*`
**nitrile** `C(#N)` => `C(=N)*`
**vinylsulfonamide** `S(=O)(=O)C=C` => `S(=O)(=O)CC*`
**bromoalkyne** `C#C[Br]` => `C(=C)*`

Additionally, `.possible_definitions` contains two (or more) that may be added experimentally (don't currently work).

**aurothiol** `S[Au]P(CC)(CC)CC` => `S[Au]*`
**aldehyde** `[C:H1]=O` => `C(O)*`

There is a quick way to get a warhead definition too

    Victor.get_warhead_definition(warhead_name)
    
In terms of Rosetta constraints (restraints), these can be added with 

    Victor.add_constraint_to_warhead(name=constrain_name, constraint=constraint)

For example:
    
* _chloroacetamide_: `AtomPair  H  145A  OY  1B HARMONIC 2.1 0.2`
* _nitrile_: `AtomPair  H  145A  NX  1B HARMONIC 2.1 0.2`
* _acrylamide_: `AtomPair  H  143A  OZ  1B HARMONIC 2.1 0.2`
* _vinylsulfonamide_: `AtomPair  H  143A  OZ1 1B HARMONIC 2.1 0.2`

Currently, only cysteine details are known to `Victor`. Cf. `.covalent_definitions`.

To convert a react_ive_ SMILES to a dummy-atom–marked react_ed_ SMILES:

    Victor.make_all_warhead_combinations(smiles, warhead_name)

### Untested backdoor
The default dummy atom can be overridden with `Monster.dummy:Chem.Mol` and `Fragmenstein.dummy_symbol:str`.
