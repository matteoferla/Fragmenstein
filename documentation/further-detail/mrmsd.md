# Combined RMSD

This submodule calculate a multiple RMSD variant which does not superpose/'align' and
bases which atoms to use on a given pairing, by default this which atoms were donors in Fragmenstein.

It requires RDKit but not pyrosetta and can be used independently of Fragmenstein core functionality.

* `mRMSD(followup: Chem.Mol, hits: Sequence[Chem.Mol], mappings: List[List[Tuple[int, int]]])` the mappings is a list of len(hits) containing lists of tuples of atom idx that go from followup to hit
* `mRMSD.from_annotated_mols(annotated_followup: Chem.Mol, hits: Sequence[Chem.Mol])` The annotated variant requires the mol to have the `_Origin` `Chem.Atom` props.
* `mRMSD.from_unannotated_mols(moved_followup: Chem.Mol, hits: Sequence[Chem.Mol], placed_followup: Chem.Mol)` the positional mapping is (re)calculated
* `mRMSD.from_other_annotated_mols(followup: Chem.Mol, hits: Sequence[Chem.Mol], annotated: Chem.Mol)` uses the second case, but `mRMSD.copy_origins(annotated, followup)` is called first.

It is a multiple RMSD, that is basically a N_atom weighted "2-mean" of RMSDs.

To properly discuss this, it is best to recap some maths.

An Euclidean distance (or 2-norm) between the vectors *a* and *b*, representing two atom positions,
is the square root of the _sum_ of the squared differences of each element/axis-position

<img src="https://render.githubusercontent.com/render/math?math=\sqrt{(a_{\overrightarrow{x}} - b_{\overrightarrow{x}})^2 + (a_{\overrightarrow{y}} - b_{\overrightarrow{y}})^2 + (a_{\overrightarrow{z}} - b_{\overrightarrow{z}})^2}">

Which can be better written

<img src="https://render.githubusercontent.com/render/math?math=\sqrt{\sum_{i=1}^{3}(a_i - b_i)^2}">

An RMSD between matrices *A* and *B*, representing two arrays atom positions, is the square root of the _average_ of the squared Euclidean distances.
If a single atom pair were compared it would nothing more than the Euclidean distance.

<img src="https://render.githubusercontent.com/render/math?math=\sqrt{\frac{\sum_{n=1}^{N} \sum_{i=1}^{3}(A_{n,i} - B_{n,i})^2}{N}}">

where _N_ is the number of atom pairs compared, _n_ is the index of a given atom pair and _i_ is simply a dimension in space (x, y, z).

So to extend the RMSD to multiple hits one can extend the pre-squared average to include all atoms pairs. One fudgey way of writing it is:

<img src="https://render.githubusercontent.com/render/math?math=\sqrt{\frac{\sum_{h=1}^{H} \sum_{n=1}^{N_h} \sum_{i=1}^{3}(A_{h,n,i} - B_{h, n, i})^2}{\sum_{h=1}^{H} N_h}}">

where _H_ is the number of hits. Note that *A* and *B* are just fudges for example purposes and they the atom pairs for one hit will differ in number between one and the next.
Written properly it would be the same as regular RMSD except A and B are the matrix concatenations for each hit

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{B} = (\mathbf{H1}_{(h\times3)} | \mathbf{H2}_{(h2\times3)} | ...)">

This means that atoms in the followup compound are re-used in the metrix as they will appear in multiple pairings.