## Overview figure notes

In the overview figure a few molecules are shown. This involved some coding as noted herein.

![overview](../../images/overview.png)

Molecules picked were random small chemicals
```python3
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

#dihydrobezene names: meta is hydroquinone, ortho is catechol, para is resorcinol
resorcinol = Chem.MolFromSmiles('c1ccc(O)cc1O')
eugenol = Chem.MolFromSmiles('Oc1ccc(cc1OC)CC=C')
catechol = Chem.MolFromSmiles('c1cccc(O)c1O')
coumarate = Chem.MolFromSmiles('c1cc(ccc1/C=C/C(=O)[O-])O')
Draw.MolsToGridImage([resorcinol, catechol, resorcinol, eugenol, coumarate])
```

Figuring out how to superpose catechol and coumarate. Rearranged index first, then just displayed atom names:

```python3
catechol = Chem.MolFromSmiles('c1(O)c(O)cccc1')
coumarate = Chem.MolFromSmiles('C(=C\\C(=O)[O-])c1ccc(O)cc1')

def mol_with_atom_index(mol):
    # https://www.rdkit.org/docs/Cookbook.html#drawing-molecules-jupyter
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

from IPython.display import display

Draw.MolsToGridImage([mol_with_atom_index( catechol ), 
                      mol_with_atom_index( coumarate )
                      ])
```

```python3
AllChem.EmbedMolecule(catechol)
AllChem.EmbedMolecule(coumarate)

manual_map = [(5, 5), (6, 0), (4, 6)]
Chem.rdMolAlign.AlignMol(catechol, coumarate, atomMap=manual_map)

import py3Dmol
view = py3Dmol.view(js="https://3dmol.org/build/3Dmol.js")
view.addModel(Chem.MolToMolBlock(coumarate), "mol", dict(style={'stick': {'colorscheme': 'cyanCarbon'}}))
view.addModel(Chem.MolToMolBlock(catechol), "mol", dict(style={'stick': {'colorscheme': 'magentaCarbon'}}))
view.zoomTo(dict(hetflag=True))
view
```

Now to verify this works with fragmenstein:

```python3
from fragmenstein import Monster
monster = Monster([catechol, coumarate]).combine()
monster.mmff_minimize()
merged :Chem.Mol = monster.positioned_mol
merged
```

Now, lets get the MCS merger

```python3
from rdkit.Chem import rdFMCS, AllChem

def superpose_by_mcs(moved: Chem.Mol, fixed: Chem.Mol) -> None:
    # superpose
    res: rdFMCS.MCSResult =rdFMCS.FindMCS([moved, fixed])
    common = Chem.MolFromSmarts(res.smartsString)
    Chem.SanitizeMol(common)
    mcs_map = list(zip(moved.GetSubstructMatch(common),
                       fixed.GetSubstructMatch(common)
                      )
                  )
    Chem.rdMolAlign.AlignMol(moved, fixed, atomMap=mcs_map)
    
superpose_by_mcs(catechol, coumarate)
monster = Monster([catechol, coumarate]).combine()
monster.mmff_minimize()
mcsed = monster.positioned_mol

Draw.MolsToGridImage(catechol, coumarate, mcsed, merged)
```

```python3
from smallworld_api import SmallWorld
display(SmallWorld.retrieve_databases())
sw = SmallWorld()
df = sw.search(Chem.MolToSmiles(merged), dist=10, db='REAL_Database_21Q3_2B')
display(df)
smilar = Chem.MolFromSmiles(df.smiles[0])
smilar
```

Save these:

```python3
for m in ('catechol', 'coumarate', 'mcsed', 'merged', 'similar'):
    Chem.MolToMolFile(globals()[m], f'test_{m}.mol')
```
The similar is horrible —the space must be hard to access.
It is also got a few tautomerisation states.

```python3
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

mol = Chem.MolFromSmiles('C=C(C)C1=CC2=C([NH]C(=O)C=C2)C(C)=C1')
enumerator = rdMolStandardize.TautomerEnumerator()
tauts = enumerator.Enumerate(mol)
Draw.MolsToGridImage(tauts)
```

In a real experiment I'd try all of them, but here there's no real protein to test.
