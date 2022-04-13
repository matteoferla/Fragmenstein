## Troubleshooting example

Say we run:
```python
from fragmenstein import Victor

victor = Victor(hits=[hit1, hit2],
      pdb_block=pdbblock)
victor.combine()
```
And it fails at the some step.

If the error is cryptic, turning on the logging to DEBUG may help.

```python
import logging

Victor.enable_stdout(level=logging.DEBUG)
```

If it fails at the `Monster` step, seeing the steps may help.
Whereas `Victor` has a `modifications` attribute,
this is filled only at the end, when `Monster` finishes.
As a result one has to look at the `modifications` attribute of `Monster` (`Dict[str, Chem.Mol]`) 
or of `Rectifier` (`List[Chem.Mol]`).

```python
from IPython.display import display

for k in victor.monster.modifications:
    print(k)
    mol = Chem.Mol(victor.monster.modifications[k])
    AllChem.Compute2DCoords(mol)
    display(mol)
```

The superscript, nominally the isotope represents the ring size.
Or seeing them in 3D:

```python
from IPython.display import display
import nglview as nv

for k in victor.monster.modifications:
    print(k)
    display(nv.show_rdkit(victor.monster.modifications[k]))
```
Note that there's no way to disable proximity bonding.
