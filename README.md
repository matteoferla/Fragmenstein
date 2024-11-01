from fragmenstein.laboratory.validator import hits_check

# Fragmenstein

## Stitched molecules
Fragmenstein: Merging, linking and placing compounds by stitching bound compounds together like a reanimated corpse.

Fragmenstein can perform two different tasks:

* **Combine** hits (merging and linking) based on their atomic overlap
* **Place** a given followup molecule based on one or more parent hits

NB. Whereas docking uses pre-generates comformers and finds the best pose that best matches the parent (if set-up to do so),
Fragmenstein creates a monstrous comformer from the parent(s) and then minimises it, optionally in the protein.
Hence why I call it a 'placement' not docking tool.

![overview](images/overview.png)

## Index

* [Installation](documentation/installation.md)
* [Theory](documentation/workings.md)
* [Command line usage](documentation/cmd.md)
* [Python usage](documentation/usage.md)
* [Examples](documentation/examples.md)
* [Limitations](documentation/limitations.md)
* [Manuscript data repository](https://github.com/matteoferla/Fragmenstein-manuscript-data)
* [Paper Authors](documentation/authors.md)
* [FAQ](documentation/FAQ.md)

## Badges and notebooks
[![Documentation Status](https://readthedocs.org/projects/fragmenstein/badge/?version=latest)](https://fragmenstein.readthedocs.io/en/latest/?badge=latest)
[![ github forks matteoferla Fragmenstein?label=Fork&style=social](https://img.shields.io/github/forks/matteoferla/Fragmenstein?label=Fork&style=social&logo=github)](https://github.com/matteoferla/Fragmenstein)
[![ github stars matteoferla Fragmenstein?style=social](https://img.shields.io/github/stars/matteoferla/Fragmenstein?style=social&logo=github)](https://github.com/matteoferla/Fragmenstein)
[![ github watchers matteoferla Fragmenstein?label=Watch&style=social](https://img.shields.io/github/watchers/matteoferla/Fragmenstein?label=Watch&style=social&logo=github)](https://github.com/matteoferla/Fragmenstein)

[![ github last-commit matteoferla Fragmenstein](https://img.shields.io/github/last-commit/matteoferla/Fragmenstein?logo=github)](https://github.com/matteoferla/Fragmenstein)
[![ github license matteoferla Fragmenstein](https://img.shields.io/github/license/matteoferla/Fragmenstein?logo=github)](https://github.com/matteoferla/Fragmenstein/raw/master/LICENCE)
[![ github release-date matteoferla Fragmenstein](https://img.shields.io/github/release-date/matteoferla/Fragmenstein?logo=github)](https://github.com/matteoferla/Fragmenstein)
[![ github commit-activity m matteoferla Fragmenstein](https://img.shields.io/github/commit-activity/m/matteoferla/Fragmenstein?logo=github)](https://github.com/matteoferla/Fragmenstein)
[![ github issues matteoferla Fragmenstein](https://img.shields.io/github/issues/matteoferla/Fragmenstein?logo=github)](https://github.com/matteoferla/Fragmenstein)
[![ github issues-closed matteoferla Fragmenstein](https://img.shields.io/github/issues-closed/matteoferla/Fragmenstein?logo=github)](https://github.com/matteoferla/Fragmenstein)

[![ pypi v fragmenstein](https://img.shields.io/pypi/v/fragmenstein?logo=python)](https://pypi.org/project/fragmenstein)
[![ pypi pyversions fragmenstein](https://img.shields.io/pypi/pyversions/fragmenstein?logo=python)](https://pypi.org/project/fragmenstein)
[![ pypi wheel fragmenstein](https://img.shields.io/pypi/wheel/fragmenstein?logo=python)](https://pypi.org/project/fragmenstein)
[![ pypi format fragmenstein](https://img.shields.io/pypi/format/fragmenstein?logo=python)](https://pypi.org/project/fragmenstein)
[![ pypi status fragmenstein](https://img.shields.io/pypi/status/fragmenstein?logo=python)](https://pypi.org/project/fragmenstein)
[![ pypi dm fragmenstein](https://img.shields.io/pypi/dm/fragmenstein?logo=python)](https://pypi.org/project/fragmenstein)

[![ codeclimate maintainability matteoferla Fragmenstein](https://img.shields.io/codeclimate/maintainability/matteoferla/Fragmenstein?logo=codeclimate)](https://codeclimate.com/github/matteoferla/Fragmenstein)
[![ codeclimate issues matteoferla Fragmenstein](https://img.shields.io/codeclimate/issues/matteoferla/Fragmenstein?logo=codeclimate)](https://codeclimate.com/github/matteoferla/Fragmenstein)
[![ codeclimate tech-debt matteoferla Fragmenstein](https://img.shields.io/codeclimate/tech-debt/matteoferla/Fragmenstein?logo=codeclimate)](https://codeclimate.com/github/matteoferla/Fragmenstein)

Example of multiple applications: [![](https://img.shields.io/youtube/views/kieDWYkzmiE)](https://www.youtube.com/watch?v=kieDWYkzmiE)

| Name                   | Colab Link                                                                                                                                                                                                     | PyRosetta | Description |
|:-----------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| :---: | :--- |
| Light                  | [![colab demo](https://img.shields.io/badge/Run_light_demo-fragmenstein.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragmenstein/blob/master/colab-notebooks/colab_playground.ipynb)  | &#10060;| Generate molecules and see how they merge<br>and how a placed compound fairs|
| Pipeline w/o Pyrosetta | [![colab demo](https://img.shields.io/badge/Run_full_demo-fragmenstein.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragmenstein/blob/master/colab-notebooks/colab_fragmenstein_Wictor.ipynb) | &#10060;| Given a template and a some hits, <br>merge them <br>and place the most similar purchasable analogues from Enamine REAL |
| Pipeline w/ PyRosetta  | [![colab demo](https://img.shields.io/badge/Run_full_demo-fragmenstein.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragmenstein/blob/master/colab-notebooks/colab_fragmenstein.ipynb) | &#10004;| Given a template and a some hits, <br>merge them <br>and place the most similar purchasable analogues from Enamine REAL |


![Ox](https://upload.wikimedia.org/wikipedia/commons/thumb/2/2f/University_of_Oxford.svg/320px-University_of_Oxford.svg.png)

## See Also

* ChemRXiv preprint — https://chemrxiv.org/engage/chemrxiv/article-details/65d751ab9138d23161b7ea38
* Fragmenstein is used in Schuller et. al. 2021
    [![SCHULLER et al](https://img.shields.io/badge/doi-10.1126%2Fsciadv.abf8711-fcb426)](https://doi.org/10.1126%2Fsciadv.abf8711)
* Figures for the upcoming manuscript are in a separate
    [repo](https://github.com/matteoferla/Fragmenstein-manuscript-data)
* The conversion of a rdkit Chem.Mol that cannot be sanitised to an analogue that can
    is done by the [molecular rectifier package](https://github.com/matteoferla/molecular_rectifier)
* The conversion of a rdkit Chem.Mol to a PyRosetta residue type (a "params file") is done via
   the [rdkit-to-params package](https://github.com/matteoferla/rdkit_to_params)
* The pipeline demo colab notebook uses Brian Shoichet's [SmallWorld webapp](https://sw.docking.org/),
    interfaced via [its API in Python](https://github.com/matteoferla/Python_SmallWorld_API)
* The playground demo colab notebook features a [JSME widget](https://github.com/matteoferla/JSME_notebook_hack) —
    [JSME](http://www.jcheminf.com/content/5/1/24) is a popular JS only molecular editor
* [Molecular Rectifier](https://github.com/matteoferla/molecular_rectifier) is used to correct the mistakes in the merged molecules,
   and is usable for other algorithms, especially de novo design via denoising diffusion probabilistic models
  (cf [blogpost discussion for the latter](https://www.blopig.com/blog/2024/09/out-of-the-box-rdkit-valid-is-an-imperfect-metric-a-review-of-the-kekulizeexception-and-nitrogen-protonation-to-correct-this/))
* Fragmenstein combine route does not check if a compound is purchasable. Above NextMove Software SmallWorld is used
  ([SmallWorld hosted by John Irwin](sw.docking.org)) to find the top N analogues,
  via [an API](https://github.com/matteoferla/Python_SmallWorld_API).
* In [Arthorian Quest](https://github.com/matteoferla/Arthorian-Quest), a parent combound is coverted with ease into an ambiguous SMARTS pattern, 
  catalogue compounds are searched with NextMove Software's Arthor ([hosted by John Irwin](arthor.docking.org)) and then placed with Fragmenstein.
* Steph Wills's [fragment network merges repo](https://github.com/stephwills/fragment_network_merges)
    enumerates superstructures of two parent hits from catalogue and places them with Fragmenstein.
* [SynDirElla](https://github.com/kate-fie/syndirella) performs a retrosynthesis of a compound, enumerates close analogues of the synthons and places their combinations
