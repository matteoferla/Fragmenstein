## Predicting new followups

Fragmenstein can also merge fragments and suggest its own 
(cf. [Moonshot result](https://discuss.postera.ai/t/fragmenstein-merging/1461)).

The method `Victor.combine` does this.

## Merging

The merging uses a 2 Å positional overlap mapping with the rings collapsed and subsequently expanded
and re-bonded by proximity with several corrections done for corner cases.

## Ring collapse

This is a seemingly "simple" solution to avoid weird bonding issues —overlaps of 5-rings with 6-rings, perpendicular rings, etc.

What happens is that all rings are replaced with a single atom that can be unpacked later.

``Ring`` class in ``core._collapse_ring`` does exactly that (inherited by ``Frankenstein``).

![collapse](../../images/atom_collapse.png)

But it can be a bit unpredictable in the upacking step after merging,
therefore it is not implemented in Victor with SMILES
—although using full merge mode is inferior to permissive unmerged mode anyway.
Instead, 

* `.collapse_ring(mol)`
* `.expand_ring(mol)`

There are two ways. remembering the bonds or making new ones by proximity for the rings.
The latter is a lot better for most cases as the rings are never 100% overlapping.

### Nitty Gritty
If stuff goes wrong and after a few checks,
one wished to debug the collapsed molecule. This is stored in `victor.fragmenstein.scaffold`:

    print(victor.fragmenstein._get_expansion_data(victor.fragmenstein.scaffold))
    
The data return is a list of length number of collapsed rings (naphthalene = 2, camphor = 2 not 3), 
here is the data for one "ringcore" atom:

    {'_ori_name': 'toluene',
     '_ori_i': -1,
     '_ori_is': '[0, 5, 4, 3, 2, 1]',
     '_neighbors': '[[1, 5], [4, 6, 0], [3, 5], [2, 4], [1, 3], [0, 2]]',
     '_xs': '[1.5, 0.75, -0.75, -1.5, -0.75, 0.75]',
     '_ys': '[0.0, 1.299, 1.299, 0.0, -1.299, -1.299]',
     '_zs': '[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]',
     '_elements': '["C", "C", "C", "C", "C", "C"]',
     '_bonds': '[["AROMATIC", "AROMATIC"], ["AROMATIC", "SINGLE", "AROMATIC"], ["AROMATIC", "AROMATIC"], ["AROMATIC", "AROMATIC"], ["AROMATIC", "AROMATIC"], ["AROMATIC", "AROMATIC"]]'}

Where `ori` is the original index. All ringcore atoms have an origin of -1.

## Logging

> See [logging_and_debugging.md](../logging_and_debugging.md) for more.

The correct logging is via Victor's journal. The logging log `Fragmenstein`.

## Valence

The class `Rectifier` attempts to fix the various issues that may have arisen.
Generally cases when the merger results in a valence that is higher than the element supports.
If so, the atom is shifted leftwards (and upwards) on the periodic table or the bond order is lowered.

## Warhead harmonisation

A issue arises merging different warheads. In which case they can be ignore, kept or the first warhead used.
Bonding to a warhead is forbidden.
Therefore, mergers may link up in unexpected ways, such as this, wherein two hits actually have different warheads.

![harmony](../../images/harmonising_warheads.png)

## Mad ones
If two rings intersect perpendicularly (_e.g._ `x0708-x2193`) the resulting bonding will be unexpected
("emergency bonding" warning appears).

![cross-ring](../../images/cross_ring.png)

I have no idea how to resolve this or whether it should be.

NB. Spiro compounds are tolerated.

## Fused rings

Even though azetine-benzene, cyclopropane-benzene and bridged compounds are chemically possible,
these are corrected as they are unlikely to be intended. As a result bridges will be removed
and for rings of size 3/4 will become 5.

## Allenes

Additionally, allene are forbidden.

## Code caveat

The code currently has many ring corrections in 'collapse ring' class (expansion part),
while bond order corrections are in 'rectifier'.
Ideally these should be merged...