## Monster entrypoints

Previously merging in `Monster` was overly complicated and an add-on, now it is cleaner as the initialisation of 
the object does not do a placement. To do that one has to do:

    monster = Monster(hits)
    monster.place(mol)
    monster.positioned_mol
    
or `place_smiles` for a SMILES.
    
While to do a merger

    monster = Monster(hits)
    monster.combine()
    monster.positioned_mol
    
Additionally, the attribute `modifications` was added/expanded, so that intermediate steps are stored
for potential inspection.
It is a dictionary, but since 3.6 dict is ordered and `.keep_copy(mol, label)` prevents over-writes.
This includes `scaffold` (the template) and `chimera` (the template with the atom symbols adapted), 
which are no longer attributes.

Equally viable alternatives are stored in a list

    monster.mol_options

`combine` still calls `simply_merge_hits`, which is used by placement too and merges by fragmentation
—unlike atom-ring absorption events. The fact that two different approaches are present is just historical.

If one wanted to merge two or more hits, independently of those in `.hits` attribute and without ring collapsing 
and rectification etc.
`simply_merge_hits` is still the method to use:

    monster = Monster([])
    monster.simply_merge_hits([molA, molB])

The `place` method accepts the argument `merging_mode`, by default it is "permissive_none",
which calls `.no_blending(broad=True)`,
but "off" (nothing), 
"full" (`.full_blending()`), 
"partial" (`.partial_blending()`)
and "none" (`.no_blending()`)
are accepted.

## Names
The names "merge" and "place" were ultra-confusingly ambiguous when it comes to the other methods, especially those of "place".

* Now all the cases where hits are partially combined for placement purposes are called "blending"
The placement of the (distorted) 3D final monster compound Victor is called "plonk".

* place and position are used synonymously
* combine and merge are used synonymously


## Future
Victor however still has the merger as a side route.

## Overlapping rings

There was a bug in the code for `Monster` that meant that if two rings from different origins that were not bonded
—i.e. without side atoms that were absorbed they were not assessed for merging. This has been corrected.

Namely, what happens to rings is controlled in `expand_ring` method, which calls `_add_novel_bonding`, where "novel"
means that does not share the same "origin" (original molecule). 
The latter method calls in turn `_get_novel_ringcore_pairs`, which now get ring marking atoms that are close or bonded.
Closeness is determined by ring atoms (not ring core markers) within 1.5 Å cutoff (viz. `_get_close_nonorigin_ringcores`)

The test case was the merger of a phenylacetic acid and phenylacetamide which overlap in two ring carbons and a terminal oxygen, resulting
in a phenylene-like ring (where one ring is an oxazine). See tests.

![phenylene](../images/phenylene.png)


    