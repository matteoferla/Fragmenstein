## Resonance

It should be said that Rosetta does not have a Drude particle,
nor would a charged virtual atom work (as they don't interact).

* If the molecule is polarisable the effects would be missed as the charges are
assigned to the topology of a molecule before placement.
* If the charge changes with isomerisation,
* either in tautomers (constitutional isomers) or stereoisomers, this will be covered by the `get_isomers` method of Victor.
* In the case of heterocycles, one may have the _possibility_ of charges running around the pi ring.

These resonance forms are not enumerated by `get_isomers`, not only because they are not
technically isomers, but because picking what a meaningful subset is tricky and
the algorithm is reported to hang in certain cases.

Indole provides a nice example of this:

```python
original = Chem.MolFromSmiles('c12c(cc[nH]2)cccc1')
#original = Chem.MolFromSmiles('c1ccccc1')

flags = Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION + \
        Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS + \
        Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS

for mol in list(AllChem.ResonanceMolSupplier(original, flags=flags)):
    AllChem.ComputeGasteigerCharges(mol)
    from rdkit.Chem.Draw import SimilarityMaps
    AllChem.ComputeGasteigerCharges(mol)
    contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
    SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
```

These are quite extreme, but they could be stabilised by bonding for example.
In fact, tryptophan has two fluorescent lifetimes, which are affected by the polarity of the solvent.

These alteranate forms are not considered by Fragmenstein.

## Pipe dream
One option to do so **would** be to score the hits
and try every resonance form until the best is found.
Then once the follow-up conformer is ready find the resonance form
with atomic charges most similar to those they map to.

This problem has two parts.
This requires the ability to change atom charge in PyRosetta
(other than remaking a pose each time, which is insane).

To do this it seems one has to generate a new residuetype,
create a residue and swap it:

So given a pose with indole (as an example):

```python
params = Params.from_smiles('c12c(cc[nH]2)cccc1')
pose = params.to_pose()
```
Get the residue type:
```python
residue = pose.residue(1)
print(residue.atomic_charge(1), residue.atom_name(1))
rt:pyrosetta.rosetta.core.chemical.ResidueType = residue.type()
```
modify it (randomly in this case):
```python
mrt = pyrosetta.rosetta.core.chemical.MutableResidueType(rt)
sac = pyrosetta.rosetta.core.chemical.SetAtomicCharge(residue.atom_name(1), 0.8)
sac.apply(mrt)
```
Get a regular residue type:
```python
nrt = pyrosetta.rosetta.core.chemical.ResidueType.make(mrt)
```
Create the replacement:
```python
new = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(nrt)
v = pyrosetta.rosetta.utility.vector1_std_pair_std_string_std_string_t(residue.natoms())
for i in range(1, residue.natoms()+1):
    v[i] = (residue.atom_name(i),residue.atom_name(i))
    new.set_xyz(i,residue.xyz(i))
```
Replace the old with the new:
```python
pose.replace_residue(1, new, False)
residue = pose.residue(1)
print(residue.atomic_charge(1), residue.atom_name(1))
```
This works, but would require correction for external bonding,
which may be a nuisance.

A partial charge is a value that goes between -1 and 1,
so finding the resonance form that best matches
would probably be a simple least squares.

The gozdillion dollar question is
how to test this is a good idea.