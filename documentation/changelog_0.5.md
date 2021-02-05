## Monster entrypoints

Previously merging in `Monster` was overly complicated and an add-on, now it is cleaner as the initialisation of 
the object does not do a placement. To do that one has to do:

    monster = Monster(hits)
    monster.place(mol)
    monster.positioned_mol
    
or `place_smiles` for a SMILES.
    
While to do a merger

    monster = Monster(hits)
    monster.merge()
    monster.positioned_mol
    
Additionally, the attribute `modifications` was added/expanded, so that intermediate steps are stored
for potential inspection.

`merge` still calls `merge_hits`.
If one wanted to merge two or more hits, independently of those in `.hits` attribute and without ring collapsing 
and rectification etc.
`merge_hits` is still the one:

    monster = Monster([])
    monster.merge_hits([molA, molB])

The `place` method accepts the argument `merging_mode`, by default it is "permissive_none",
which calls `.no_merging(broad=True)`,
but "off" (nothing), 
"full" (`.full_merging()`), 
"partial" (`.partial_merging()`)
and "none" (`.no_merging()`)
are accepted.

## Future
Victor however still has the merger as a side route.

The names "merge" and "place" are ambiguous when it comes to the other methods, especially those of "place".
"position" is also used for the (distored) 3D final compound.

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


    