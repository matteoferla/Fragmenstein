## Three modes

As [discussed](README.md), it has three main parts:
* the pipeline (Victor class)
* the minimiser (Egor class) and
* the followup placing algorithm (Fragmenstein class)

The latter has three modes as discussed in [the fragmenstein page](fragmenstein.md): "full-merging", "partial-merging" and "no-merging" for mapping.

Each has advantages and disadvantage. For the [covid](covid.md) data, an analysis was done to figure out which.

50 followups were randomly chosen and placed with the three modes (165 total poses were done, +10% for safety).

The medians of the various scores are:

| mode    |      ∆∆G [kcal/mol] |   concatenated RMSD [Å] |   Number of constrained atoms |   Runtime [s] |
|:--------|---------:|----------:|----------------------:|----------:|
| full    | -3.55336 |  0.940974 |                    12 |   16.6215 |
| none    | -4.02654 |  1.0435   |                    14 |   17.5834 |
| partial | -3.23976 |  0.970151 |                    13 |   16.7117 |

No merging mapping is better for both number of atoms used and the ∆∆G_bind.
It scored worse for time and slightly for RMSD.

![compare-N-atoms.png](images/compare-N-atoms.png)
![compare-rmsd.png](images/compare-rmsd.png)
![compare-∆∆G.png](images/compare-∆∆G.png)
![compare-runtime.png](images/compare-runtime.png)
It is clear that the no merging mapping when it fails it fails badly.
So let's look solely at all ∆∆G < 0 ones:
![compare-∆∆G2.png](images/compare-∆∆G2.png)

The outliers in ∆∆G are caused by erroneous hits being claimed to be inspirations.
Here is one of these, WAR-XCH-b6889685-63 in the same colours as the graphs:

![image](images/WAR-XCH-b6889685-63.png)

Namely the outlier is chosen as an inspiration and the bond length is not violated as the atoms from the mapped parts are 
not neighbours, so the 3 Å limit does not apply. Curiously, full-merger mapped to one cluster, while partial merger to the loner.

The time problem is substantial for the none.
![compare-timeproblem.png](compare-timeproblem.png).

These are the offenders:
|     | name                | mode    |       ∆∆G |   comRMSD |   N_constrained_atoms |   runtime |   N_hits |
|----:|:--------------------|:--------|----------:|----------:|----------------------:|----------:|---------:|
|  36 | JAG-SYN-9c2cd0bd-13 | none    | nan       | nan       |                    18 | 2164.83   |        6 |
| 156 | WAR-XCH-b6889685-20 | none    |  -5.81899 |   1.43685 |                    10 |  528.47   |        5 |
|  93 | CHR-SOS-1f323c23-1  | none    |  -7.99253 |   0.65989 |                    16 |   85.4308 |        5 |
|  98 | RAM-SYN-2a37ce6c-4  | partial |  -8.13456 |   2.2751  |                    11 |   79.835  |        4 |
| 153 | DRA-CSI-47e38074-1  | none    |  -4.02654 |   1.0435  |                    13 |   67.3437 |        6 |

The `JAG-SYN-9c2cd0bd-13` has nothing wrong with it. It dies after the igor step due to sanitisation (fixed).
The "none" mode excluded the following `['x0689', 'x0691', 'x0692', 'x0705', 'x0708']`,
As a human at first glance I would have excluded `['x0689', 'x0692', 'x0705', 'x0708']`, so keeping `['x0691', 'x0731']`,
but looking better at it, `'x0691'` would be in the list. So perfect!

However, 33 minutes is not an acceptable runtime.
