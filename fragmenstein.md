# Fragmenstein
## Description

The followup placing algorithm (`Fragmenstein` class) has three modes.

The three modes rely heavily on mapping one-to-one atomic coordinates of overlapping atoms between hits
within a 2Å radius (see (Position Mapping code)[fragmenstein/core/positional_mapping.py]).

The three have different strengths. And which is better depends actually on the dataset.

The argument `merging_mode` (`full`|`partial`|`none`|`off`) chooses which is used.

## Full merging

_The main documentation for this mode can be found (here)[fragmenstein_full.md]_.

The inspiration hits are merged creating a new merged template and the followup is mapped onto this.
This makes it much faster than the unmerged mode, so sliding scale of MCS mappings is done
—a very strict MCS mapping is done, then a series of MCS ranging from very lax to more strict are done until one is found that supersets the strictest one.

This was the primary method used, until partial mapping was introduced to overcome the incorrect hit issue.

### Pros

* Obedient
* Can be used to suggest new followups

### Pro/Con
* Code to stop mapping of unconnected hits (`fuse`) disabled/buggy, because of erroneous hits are a bit issue.

### Cons

* Slower than partial mapping (reason unclear)
* Sensitive to too many and to incorrect hits (a COVID moonshot dataset problem)
* Ring overlaps can result in really odd rings(&lowast;)

&lowast; See [work in progress](wip.md) for the ring collapsing code which aims to fix this.

## Partial mapping
As above but inconsistent hits are excluded.
The "inconsistent" hits (`dodgies` in the code) are those who when seen in trio of hits do not map consistently.
This speeds up the code and avoids the pitfalls of incorrect hits and too many hits.

### Pros
* Faster than full
* Less prone to too many hits/wrong hits

### Cons
* often behaves differently than intended by ignoring good hits

# No merging
This method (self standing in `core.unmerge_mapper.py`) tries to maps each hit to the compound.
The maps must satisfy some conditions, two mapped atoms cannot occupy the same 2 Å overlapping position
and no bond in the map can be over 3 Å long —if the atoms in between are not placed there is no constraint (which is a problem).
This mills through all possible MCS pairing of the various inspiration hits, so the time taken increases to the power of the hits.

### Pros
* Better at mapping

### Cons
* Slower
* Vulnerable to erroneous hits

## Comparison
For a comparison see [three modes compared](three_modes_compared.md).

