from fragmenstein import Laboratory

# Pythonic usage

### Prerequisites

**The parent hits have to be 3D embedded in the same reference frame**,
that is, they have to be extracted from superposed crystal structures or docked in the same receptor.
If the parent hits are not in the same reference frame, the results will be rubbish and fail
with either distance error or fail to minimise as the monster is jammed in the middle of a protein chain.

The receptor/template protein has to lack a ligand in the pocket of interest unless this is expected to be there.
For example, if the protein has ATP in the chosen reference structure,
but the structure of the parent hit lacks this because it is displaced, then it needs to be absent
—it will be automatically stripped unless an ion or a canonical amino acid, but is best to not take chances.
If the ligand should be there, for example a cofactor like PLP, then it should be present in the parent hit,
and the topology file provided (vide infra).

### Classes

There are four main classes —named after characters from the Fragmenstein book and movies:

* `Monster` makes the stitched together molecules indepently of the protein — [documentation](monster/monster.md)
* `Igor` uses PyRosetta to minimise in the protein the fragmenstein monster followup — [documentation](further-detail/igor.md)
* `Victor` is a pipeline that calls the parts, with several features, such as warhead switching —[documentation](further-detail/victor.md)
* `Laboratory` does all the combinatorial operations with Victor (specific case)

NB. In the absence of `pyrosetta` (which requires an academic licence), all bar ``Igor`` work and 
alternative Victor classes need to be used, for example 
`Wictor` (RDkit minimisation only), `OpenVictor (using OpenMM).

Additionally, there are a few minor classes.

One of these is ``mRMSD``, a multiple RMSD variant which does not superpose/align and bases which atoms 
to use on coordinates —[documentation](further-detail/mrmsd.md)

The class `Walton` performs geometric manipulations of compounds, to set them up to demonstrate
features of Fragmenstein (like captain Walton, it does not partake in the plot, but is key to the narration)

There are two module hosted elsewhere:

* ``Rectifier`` from [molecular_rectifier](https://github.com/matteoferla/molecular_rectifier) is a class that corrects mistakes in the molecule automatically merged by ``Monster``.
* ``Params`` from [rdkit to params module](https://github.com/matteoferla/rdkit_to_params) parameterises the ligands

### Combine
It can also merge and link fragment hits by itself and find the best scoring mergers.
For details about linking see [linking notes](further-detail/linking.md).
It uses the same overlapping position clustering, but also has a decent amount of impossible/uncommon chemistry prevention.

Monster:

```python
from fragmenstein import Monster
hits: List[Chem.Mol] = ... # 1 or more RDKit.Chem.Mol, sanitised, w/ conformer, preferably without explicit Hs
monster = Monster(hits=hits)
monster.combine()
monster.positioned_mol #: RDKit.Chem.Mol
```

Victor:

```python
from fragmenstein import Victor
import pyrosetta
pyrosetta.init( extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

hits: List[Chem.Mol] = ...
victor = Victor(hits=hits, # List of 1 or more RDKit.Chem.Mol
            pdb_filename='foo.pdb',  # or pdb_block='ATOM 1 MET ...'
            covalent_resi=1) # if not covalent, just put the first residue or something.
victor.combine()
victor.minimized_mol
```
The PyRosetta init step can be done with the helper function:
```python
Igor.init_pyrosetta()
```

The two seem similar, but Victor places with Monster and minimises with Igor.
As a result it has ∆G_bind energy score (difference between holo minus apo+ligand Gibbs free energy predictions)
(`victor.ddG`).

Fragmenstein is not really a docking algorithm as it does not find the pose with the **lowest energy** 
within a given volume.
Consequently, it is a method to find how **faithful** is a given followup to the hits provided.
Hence the minimised pose should be assessed by the RMSD metric or similar
and the ∆G_bind score used solely as a cutoff —lower than zero.

For a large number of combination:

```python
from fragmenstein import Laboratory

pdbblock: str = ... # a PDB block
hits: List[Chem.Mol] = ... # 1 or more RDKit.Chem.Mol
lab = Laboratory(pdbblock=pdbblock, covalent_resi=None)
combinations:pd.DataFrame = lab.combine(hits, n_cores=28)
```

### Place
Here is [an interactive example of placed molecules](https://michelanglo.sgc.ox.ac.uk/r/fragmenstein).

It is rather tolerant to erroneous/excessive submissions (by automatically excluding them)
and can energy minimise strained conformations.
![summary](../images/new_summary.jpg)

Three mapping approaches were tested, but the key is that hits are pairwise mapped to each other by means 
of one-to-one atom matching based upon position as opposed to similarity which is easily led astray. 
For example, note here that the benzene and the pyridine rings overlap, not the two pyridine rings:

<img src="../images/position_over_mcs.jpg" width="300px">

### RDkit only and OpenMM

PyRosetta is needed for the pocket-centric minimisation.
Two alternatives are available:

* `Wictor` (without): stops at the RDKit minimisation
* `OpenVictor` (with OpenMM): uses OpenMM to minimise in the protein

Whereas the PyRosetta steps operate via Igor, OpenVictor uses Fritz.
OpenMM is a lot slower than PyRosetta on CPU only,
but is free, open source and potentially more accurate.

Igor is a much larger class as it needs to disable rotamer sampling and other things,
which is not an issue in OpenMM.

A further detail is that openMM is already parallel,
therefore when using with `Laboratory` request only one core.
```python
from fragmenstein import Laboratory, OpenVictor
Laboratory.Victor = OpenVictor

lab = Laboratory(pdbblock=MPro.get_template())
combinations: pd.DataFrame = lab.combine(hits,
                                     n_cores=1,  # 1 core unless $OPENMM_CPU_THREADS is set
                                     timeout=600,  # 2 minutes
                                     combination_size=2,  # pairwise
                                     max_tasks=0)  # 0 is no chunking
```

## Score

The ∆Gbind and RMSD do go in tandem —if something is bad it is bad for both—,
but for ranking of valid compounds more is needed.
Firstly, ∆Gbind is misled by high number of unconserved (novel) atoms.
Secondly, it does not reflect how many interactions are kept.
Thirdly, it favours larger compounds (cf. ligand efficiency discussion).

Therefore, there is a multiparametric penalty score (negative = good) calculated based on
arbitrary weights and called `ad_hoc_penalty` in the output pd.DataFrame.

```python
import pandas as pd
from fragmenstein import Laboratory, unbinarize, cli_default_settings
import warnings, operator
from rdkit import rdBase

# ## Weights
# let's arbitrarily set the weights
weights = cli_default_settings['weights']
# thi is a probability scaled sum of the interactions: unique interactions will give bigger numbers
weights['interaction_uniqueness_metric'] = -0. # formerly -5
weights['N_interactions_lost'] = 5
weights['N_unconstrained_atoms'] = 0.5

# ## DataFrames
# placement results of the frag-hits 'redocked':
hit_replacements: pd.DataFrame = ...  
# let's pretend someone one has lost the list of hits, not ideal, but it happens
# getting the after going through the cycle
hits = hit_replacements.hit_binaries.apply(operator.itemgetter(0)).apply(unbinarize).to_list()
df: pd.DataFrame = ... # results

# ## Score
with rdBase.BlockLogs():  # suppresses RDKit warnings
    Laboratory.score(df, hit_replacements,
                     suffix='_manual', hits=hits, weights=weights
                     )

df = df.sort_values('ad_hoc_penalty').reset_index(drop=True).copy()
```

The `cluster` column is a Butina clustering of Tanimoto distances by Morgan fingerprints (cutoff=0.35 by default).

If the frag-hits were spread out, no kinetic data is available, the target is a protein-protein interaction,
then clustering by interactions and picking the best within each is a good idea. See below.

## Interactions

Victor uses PLIP to predict interactions between the protein and the ligand.
(See installation for note on issues with this).

```python
vicky = Victor(hits=hits, pdb_filename='receptor.pdb')
itxns: Dict[Tuple[str, str, int]] = victor.get_plip_interactions()
```
In Laboratory wrapper:

```python
lab = Laboratory(pdbblock=pdbblock, run_plip=True)
```

The key is a tuple of interaction type (str), residue type (str) and PDB residue index (int).
This can be a bit ugly when exporting to CSV:

```python
import pandas as pd

placed: pd.DataFrame = ...
joined_cols = {col: ':'.join(map(str, col)) for col in placed.columns.to_list() if
               isinstance(col, tuple) and placed[col].sum() > 0}
other_keepers = ['name', 'smiles', 'ad_hoc_penalty', '∆∆G', 'comRMSD', 'hit_names', 'cluster', 'N_interactions',
                'N_constrained_atoms', 'N_unconstrained_atoms', 'N_interactions_kept', 'N_interactions_lost',
                'interaction_uniqueness_metric', 'max_hit_Tanimoto', 'N_PAINS', 'strain_per_HA', 'MW', 'HAC',
                'percent_hybrid', 'largest_ring', 'N_rotatable_bonds']
placed.loc[(placed.ad_hoc_penalty < 0) & (placed['∆∆G'] < -5) & (placed.MW > 250)] \
        .sort_values('ad_hoc_penalty') \
        [other_keepers + list(joined_cols)] \
        .rename(columns={'∆∆G': 'dG_bind', 'ad_hoc_penalty': 'combined_score', **joined_cols})
        .to_csv('fragmenstein_results.csv')
```

To cluster by interaction, play around with the method as there are several options,
which may be dependent on the data, which is zero inflated so biases abound and 
something that should be pretended to be simple will have a lot of complexity brushed under the carpet.
Here is an example of two different approaches.

The interactions are actually not really boolean, as a virtual compound may form more than one interaction,
but this is uncommon. The interactions are also highly covariant.
The types of interactions are not equal: a hydrophobic is worth less energetic than a salt bridge.

So one can treat the data as binary (KModes) or use KMeans and pretend there is no covariance and zero inflation.
Here is a simple multiple option snippet:

```python
import pandas as pd
from scipy.cluster.vq import kmeans, vq
from kmodes.kmodes import KModes
import pandera.typing as pdt
from collections import defaultdict
from Bio.SeqUtils import seq1
import enum


class clustering_algo(enum.Enum):
    KMEAN = 1
    KMODE = 2


class scaling_used(enum.Enum):
    NONE = 0
    CRUDE = 1
    PROB = 2
    FF = 3


# rubbish scaled
def rubbish_scale(col):
    """
    This is in effect inferior intermolecular LJ
    A hydrophobic is worth half, a salt bridge 2x
    
    :param col: 
    :return: 
    """
    if col.name[0] == 'hydroph_interaction':
        w = 0.5
    if col.name[0] in 'halogenbond':
        w = 1.5
    elif col.name[0] == 'saltbridge':
        w = 2.0
    else:
        w = 1.0
    return col * w

class Ranker:
    
    def __init__(self):
        self.rank = defaultdict(int)
        
    def __call__(self, c):
        self.rank[c] += 1
        return self.rank[c]

def munge_data_for_clustering(analogs):
    intxn_cols = [c for c in analogs.columns if isinstance(c, tuple)]  # assumes itxn cols are tuples
    data_for_clustering = analogs.loc[analogs.outcome == 'acceptable'][intxn_cols].fillna(0).copy()
    return data_for_clustering

def calc_labels(data_for_clustering, scaling_mode):
    if scaling_mode == scaling_used.PROB:
        # probability scaled
        tallies = data_for_clustering.sum().to_dict()
        data_for_clustering = data_for_clustering.apply(lambda col: col / tallies[col.name], axis=0).fillna(0)
        algo = clustering_algo.KMEAN
    elif scaling_mode == scaling_used.CRUDE:
        # A hydrophobic is worth half, a salt bridge 2x
        data_for_clustering = data_for_clustering.apply(rubbish_scale, axis=0).fillna(0)
        algo = clustering_algo.KMEAN
    elif scaling_mode == scaling_used.FF:
        raise NotImplemented
    elif scaling_mode == scaling_used.NONE:
        # no scaling so Kmode
        algo = clustering_algo.KMODE
    else:
        raise ValueError
    
    if algo == clustering_algo.KMEAN:
        centroid, variance = kmeans(data_for_clustering.values, k)
        labels, _ = vq(data_for_clustering.values, centroid)
    elif algo == clustering_algo.KMODE:
        km = KModes(n_clusters=k, init='Huang', n_init=5, verbose=1)
        labels = km.fit_predict(data_for_clustering)
    else:
        raise ValueError
    return labels

# ---------------
k = 10  # number of clusters sought
scaling_mode = scaling_used.PROB
analogs: pd.DataFrame = ...
data_for_clustering = munge_data_for_clustering(analogs)
labels = calc_labels(data_for_clustering, scaling_mode)
# list to series first for the correct indices:
data_for_clustering['cluster_label']: pdt.Series[int] = labels
analogs['cluster_label']: pdt.Series[float] = data_for_clustering.cluster_label
analogs['cluster_label']: pdt.Series[int] = analogs['cluster_label'].fillna(-1).astype(int)
analogs = analogs.sort_values('full_penalty').drop_duplicates('cluster').reset_index(drop=True)
analogs['cluster_rank'] = analogs['cluster_label'].apply(Ranker())
analogs = analogs.sort_values(['cluster_rank', 'full_penalty']).reset_index(drop=True).copy()
```

## E-amide isomerism

The static method `Monster.inspect_amide_torsions` called on a Chem.Mol will return a list of the amide torsions.
This is because some mergers do result in cis (E) amides, which are not ideal.
By default, the MMFF forcefield will penalise exocyclic cis-amides.
There is a setting, `ff_prevent_cis` (as env `$FRAGMENSTEIN_FF_PREVENT_CIS`),
which when set to `false` will disable the addition of a penalty against cis-amides to the MMFF forcefield
when `Monster.mmff_minimize` is called.