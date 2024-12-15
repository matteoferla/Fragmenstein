# Installation

* Essential: Fragmenstein
* Optional: Pyrosetta or OpenMM, PLIP and OpenBabel

### Fragmenstein and dependencies

Python 3.6 or above. It can be install from pipy

    python -m pip install fragmenstein

### Requires Pyrosetta

PyRosetta is optional, but it is required for full operations.
Without it the pocket sidechains are not minimised, which requires a perfect receptor/template/protein.

:warning: PyRosetta no longer runs on CentOS 7 due to old kernel headers (cf. [blog post](https://blog.matteoferla.com/2022/11/glibc-236-vs-centos-7-tale-of-failure.html)).

:warning: PyRosetta is not available on Windows.
For that you need to [install Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
and install PyRosetta on the Linux subsystem.

#### Modern way

This no longer requires a password:

```bash
pip install pyrosetta-installer 
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
```

#### Details from former way

Pyrosetta used to require a password to be downloaded (academic licence) obtained by https://els2.comotion.uw.edu/product/pyrosetta. 
This is a different licence from the Rosetta one. The username of the Rosetta binaries is formatted variant of "academic user", 
while the PyRosetta is the name of a researcher whose name bares an important concept in protein folding,
like boltzmann + constant (but is not that). 
Pyrosetta can be downloaded via a browser from http://www.pyrosetta.org/dow. Or in the terminal via:

    
```bash
curl -u 👾👾👾:👾👾👾https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python38.linux/PyRosetta4.Release.python38.linux.release-NNN.tar.bz2 -o a.tar.bz2
tar -xf a.tar.bz2
cd PyRosetta4.Release.python38.linux
sudo pip3 install .
```

or using conda

or using `install_pyrosetta` from the `pyrosetta-help` package.

```bash
pip install pyrosetta-help
PYROSETTA_USERNAME=👾👾👾 PYROSETTA_PASSWORD=👾👾👾 install_pyrosetta
```
The `PYROSETTA_USERNAME` and `PYROSETTA_PASSWORD` are environment variables,
which should not be shared publicly (i.e. store them as private environmental variables
in your target application).

## PLIP

PLIP is used for the interactions and it is optional.

It is involved by Victor method `get_plip_interactions` and Laboratory will use it if it finds it.

To install PLIP, sometimes issues happen as openbabel is installed via conda, so I do this:

```bash
conda install -c conda-forge openbabel lxml
conda install -c conda-forge --no-deps plip
```

