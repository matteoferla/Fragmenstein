## History

The code was originally written to make a conformer by merging as opposed to by mapping each combination.

The automerging route was therefore built on top of it.

The `_VictorAutomergeMixin` contains the automerging code,
most notably the classmethod `combine`.
This uses `merge_hits` from the main class, written first.

Then to address the problem of rings merging oddly, the ring collapse functionality was added
and can be found in `_collapse_ring`.

Then the ring rejoining algorithm was tweaked to bond by default by proximity and not memory.
Now the rings are bonded by memory if original, by proximity if novel.

Then to address odd bonding, the class `Rectifier` was added.
This is a class called by Fragmenstein and not inherited.

However, the ring collapse code deals with odd ring sizes.

## Mapping
Also, mapping is done by the class `GPM` in `positional_mapping.py`.
This is used by both SMILES no-merging mode and the automerging mode.