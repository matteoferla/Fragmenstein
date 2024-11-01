## Number of parent hits

The algorithm is not flexible in respect to parent hits: it does not arbitrarily mix and match:
the user needs to do that.
There are some steps taken to discard atoms that are red herrings etc., but these are limited.
If you have a compound and don't know which of a 100 parent hits was the parent hit,
place iteratively the compound and pick the best, don't run all 100 parent hits as a single job.

## Catalogue and fragment sociability

Whether a compound has many superstructures on one or more expansion vectors
is called fragment sociability (cf. Astex paper).
Fragmenstein is heavily affected by this as it is a generative method.

Chances are, for a given combined compound, even if it looks totally obvious to make from common building blocks,
it is unlikely to be found in catalogue space.
The bigger the compound is, the more likely it is that it will require a custom synthesis.

Often it is rather non-sensical. For example, a certain vendor that starts with an E does not have many 1,4-disubstituted triazoles (copper AAC) and very few 1,5 ones (ruthenium AAC),
yet the reaction is so common and facile that Huisgen won the 2016 Nobel prize for it along with Heck and Suzuki.
Another example, the filters in place might preclude some groups like carboxylic acids.

This is why a search step is needed.

Alternatively, I used [Arthor](arthor.docking.org) via [Arthorian-Quest](https://github.com/matteoferla/Arthorian-Quest),
by making an ambiguous SMARTS pattern and explore what analogues work,
or see what Reaxys suggests I do and try and fix the nasty bits.

## Angles

A problem that sometimes arises is that the angles are not right,
when linking two arenes.

* Do they need to be aromatic? Or would a greasy 3D compound work, like norbornane?
* Can one be picked to be reduced and the other kept aromatic?
* Why not look for a shape and colour analogue...







