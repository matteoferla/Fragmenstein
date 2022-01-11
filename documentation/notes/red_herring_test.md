## Red Herring test

One of the tests fails to exclude the red herring.

For MPro 2_ACL the fragments stated where x1249, x0305 and x0692.
But x1249 is the red herring. However, Fragmenstein falls for it.

To display the hits, using code copied 
from [Matteo's blog](https://blog.matteoferla.com/2021/02/multiple-poses-in-nglview.html):

```python3

import nglview as nv
from rdkit import Chem
from io import StringIO
from typing import *
from warnings import warn

view = nv.NGLWidget()

def get_ggplot_colour_scale(n:int=7) -> Iterable[NewType('ColorHex', str)]:
    ggplot_color_scales = {1: ['#F8766D'],
                           2: ['#F8766D', '#00B4C4'],
                           3: ['#F8766D', '#00BA38', '#619CFF'],
                           4: ['#F8766D', '#7CAE00', '#00BFC4', '#C77CFF'],
                           7: ['#F8766D', '#C49A00','#53B400','#00C094','#00B6EB','#A58AFF','#FB61D7']
                           }
    if n in ggplot_color_scales:
        return iter(ggplot_color_scales[n])
    else:
        return iter(ggplot_color_scales[7])
        

def crete_multiple_view(*mols: Chem.Mol) -> nv.NGLWidget:
    if len(mols) == 1 and isinstance(mols[0], Sequence):
        warn('Expected each `Chem.Mol` as an argument, got  `Sequence[Chem.Mol]`')
        mols = mols[0]
    colors = get_ggplot_colour_scale(len(mols))
    view = nv.NGLWidget()
    for m, mol in enumerate(mols): #: Tuple[int, Chem.Mol]
        fh = StringIO(Chem.MolToMolBlock(mol))
        view.add_component(fh, ext='mol')
        view.clear_representations(component=m)
        view.add_licorice(colorValue=next(colors), component=m, multipleBond='symmetric')
    return view

def show_colors(*hit_names: str) -> None:
    if len(hit_names) == 1 and isinstance(hit_names[0], Sequence):
        warn('Expected each `str` as an argument, got  `Sequence[str]`')
        hit_names = hit_names[0]
    from IPython.display import display, HTML
    template = '<span style="color: {1}">{0}</span>'
    zipped = zip(hit_names, get_ggplot_colour_scale(len(hit_names)))
    spans = [template.format(name, color) for name, color in zipped]
    display(HTML('<br/>'.join(spans)))
```

<span style="color: #F8766D">x1249</span><br/>
<span style="color: #00BA38">x0692</span><br/>
<span style="color: #619CFF">x0305</span>

![red](../../images/red-herring.png)

Molecule to make:

```python3
Chem.MolFromSmiles('CCNc1ncc(C#N)cc1CN1CCN(C(=O)C*)CC1')
```

The hit x0305 is a basically a substracture, but the ring on x0692, clashes with it
as the latter's methyl substituent in para and the nitrile substituent in para of the former
compete for the same mapped atom, while x1249 does not, so scores better, even if it matches less.

The automatic elimitation machinery is cool for mix and matching, but it is poor at predicting human
intentions, which means it is not actually a substitute.

