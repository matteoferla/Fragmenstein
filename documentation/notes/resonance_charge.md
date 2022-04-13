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