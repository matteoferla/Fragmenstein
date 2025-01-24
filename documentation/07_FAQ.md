## It is failing

* If you get a load of cryptic warnings: this is normal if you did not configure the logger (`Victor.journal`). Make sure to set `Victor.enable_stdout(logging.WARNING)` or at the very least call `Victor.capture_rdkit_log()`
* Your hits are not in the same reference frame as your protein or each other
* Your hits must be sanitised (`AllChem.SanitizeMol`)
* Your hits cannot be None
* If there are two or more disconnected molecules in your hit they will be joined. So be sure you do not have 2+ copies of your hit in you Chem.Mol object (e.g. extracted from a homodimeric structure)?
* The hits have to have unique titles ('_Name' property in RDKit)
* The hits are best without explicit hydrogens. Unlikely to cause issues. The hydrogens are not removed because there is a case when it is needed.
* Your hits have to have bond order (use Openbabel to assign it) otherwise you will get alicyclics

## I do not know what parent hits to use

Do not provide all your hits for either combine or place as it will take forever.

Fragmenstein can handle more than 2 parent fragments in both routes: this is done via iterative operations on pairs of hits. 
Due to the exhaustive matching approach, inference time increases exponentially, limiting the merging of all hits into a single compound or placing a compound against all fragment hits. 
Depending on the user’s choice, it either tries to use all fragments, raising an error when a fragment does not “fit,” or disregards fragments.
Combining multiple fragments. In the combination route, it first identifies, in a pairwise fashion, the overlapping atoms of the hits, then it mergers the pair of hits with the most overlaps and then merges the merged compound with the next hit and so forth. The order in which the parent compounds were provided matters as it influences the priority in the case of ties, and, among other things, what element is chosen: for example the merger of an superposed furan and an unsubstituted benzene will yield a furan or a benzene depending on the order (cf. Figure 1B). However, a different approach is preferable regardless. For mergers/linkers, catalogue space (for either final compounds or building blocks) can be limiting, especially with compounds having three or more fused rings, so a tournament approach, going through catalogue space, is advisable.
Placing multiple fragments. In placements, the iteration processed by identifying first the parent compound with the largest shared substructure, and then iteratively finding the next compound with the most further substructures that are consistent with the already mapped structure and are novel additions. Where the design history is unfortunately lost, iterating pairs and assessing the results is more effective, because in cases where three or more parent hits are used it is often the case where several parents contribute a single substituent atom, whereas a single or a pair of hits contribute the scaffold. A special case is present in scaffold hopping, where the user removes the core from the parent hit and places a compound with a different scaffold. Therein ambiguous single atoms like multiple methyl groups are best removed by the user, but will likely be flagged as non-productive and ignored, effectively making it no more than a three-way placement.

## I want to do something complicated

There are a few shortcuts I can recommend.
Although, I should say that if the compound you are trying to make is complicated, then it most likely will not work IRL.

It is possible to specify which atoms to map (`custom_map`) with `Monster`, `Victor` and `Laboratory`.
But this can be complicated as it requires figuring out which atoms are which.
A shortcut is to use `AllChem.DeleteSubstructs` to remove the parts from the parent hits.
This can be used nicely with `AllChem.GetMolFrags`.
Say the parent hit has an amide and I want the part which includes the zeroth atom:

    fore: Chem.Mol = AllChem.GetMolFrags(AllChem.DeleteSubstructs(hit, Chem.MolFromSmiles('O=CN')), asMols=True)[0]


Whereras you can use monster recursively:

    intermediate: Chem.Mol = Monster([hit1, hit2]).combine().positioned_mol
    
    mod = Chem.RWMol(AllChem.RemoveAllHs(intermediate))
    mod.RemoveBond(43, 65)
    mod.AddBond(0, 65, Chem.BondType.SINGLE)
    neo: Chem.Mol = mod.GetMol()

    final: Chem.Mol = Monster([neo, hit3]).combine().positioned_mol    

You might get weird errors, especially if using `Chem.RWMol` multiple times.
If you are getting something weird do a round trip:

    good_mol: Chem.Mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(naughty_mol))

## Linkers

> How are linker atoms chosen?  How does the chemical rectification work?  

This is covered in detail in [linking.md](further-detail/linking.md)



## Wrong tautomer

If a tautomeric molecule were placed with Fragmenstein in a given tautomer,
but the binding mode best suited the other tautomer and the latter mimicked the parent hit,
the binding mode of the latter ought to be seen, but with the former tautomer.
It is pure conjecture as no such test cases were found in a manual search, 
whereas the unfortunately common scenario where a user docks a carboxylic acid, yet the bound form is a carboxylate, 
is circumvented in Fragmenstein if the parent hit behaved like the carboxylate,
but ought to have been circumvented by better ligand prep.

NB. `Laboratory().place(..., expand_isomers=True)` can expand stereoisomers, but does not expand tautomers.

## Fragmenstein and sterically hindrance

> How does Fragmenstein avoid sterically unfavourable bond vectors, such as axial CH in cyclohexane or the conserved amide NH in the Covid Moonshot isoquinolinyl amides, for structural elaboration. 

In Fragmenstein the rules governing mergers prevent impossible chemistry and only marginally unfeasible chemistry, 
majorly excess of 3-membered rings. 
However, in the case of an atomic overlap between atoms of bulky substituents one might be ‘absorbed’ into the other resulting in either a bridged alicylic or a single substituent as there is a rule to prevent propallanes. 
In terms of forbidden torsions, an explicitly encoded case are exocyclic secondary cis-amides: these are prevented during minimisation: were an aniline merged to a benzoic acid to form (E)-phenylbenzamide it would be forced into Z if possible or rejected. 

If a user was intent on preserving a substructure, like an unsubstituted amine, they could flag the atoms with the ‘protection’ mechanism: this is used automatically for warhead/reaction-product moieties, 
but is an RDKit property on the `Chem.Mol` object of the parent hit.

## Large binding sites

> How does the code behave for large binding sites?

Fragmenstein fails more potential pairings or placed compounds when the cavity is restrictive, a common scenario is two adjacent perpendicular arenes. With the caveat that such cases are more suitable for scaffold hopping than merging in the first place, Fragmenstein will likely try to dearomatise one ring or both, yielding a fused aromatic-alicyclic or a full alicyclic, which generally fail the RMSD filter.
When the region of interest is large, as is often the case for protein-protein interfaces, in a benchmarking scenario, there will be an overwhelming number of acceptable purchasable analogues, which appears as a positive. In an applied scenario, the fragment-hits are often many, but distributed heterogeneously and weak-binding: certain hypotheses/series will need to be addressed by complementary methods, such as catalogue enumeration of superstructures to join two moieties over 5Å apart, which can however be placed with Fragmenstein.

