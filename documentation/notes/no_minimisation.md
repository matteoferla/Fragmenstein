## No PyRosetta

> This note will quickly become outdated once Fritz enters the scene.

Fragmenstein can import without PyRosetta, but it will not be able to perform any minimisation.
Only `Monster` will work.
If `Igor` is imported, it will not be able to perform any minimisation 
and crash with the final outputs being a bit hidden.

```python
from fragmenstein import Victor  # RuntimeWarning: PyRosetta is not installed. A mock object is loaded. Any calls will fail.
from fragmenstein import __version__ as fragn_version

print(fragn_version)  # 0.9.12.6
```

A way around it is to create a subclass of `Victor` that overrides the methods that call `Igor`.

```python
from rdkit import Chem

class Wictor(Victor):
    """
    This Victor does not call Igor
    """
    def _calculate_combination_thermo(self):
        # override igor.
        pass
    
    
    def _calculate_placement_thermo(self):
        # override igor.
        pass
    
    def post_monster_step(self):
        # this is a black overridable methods that will be the last thing called
        self.minimized_mol: Chem.Mol = self.monster.positioned_mol
        self.minimized_pdbblock: str = self.unminimized_pdbblock
```

```python
Wictor = ... # see above

from fragmenstein.demo import Mac1

Wictor.capture_rdkit_log()
#Wictor.enable_stdout()
Wictor.error_to_catch = ()
wicky = Wictor(hits=Mac1.get_n_filtered_mols(2),
               pdb_block=Mac1.get_template(),
              )

wicky.combine()
wicky.to_nglview()
```


