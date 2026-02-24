"""
See documentation/notes/colours.md

``divergent_colors`` are R-style CIE HChL divergent.

These were generated via:

.. code-block:: python
     !pip install git+https://github.com/retostauffer/python-colorspace
     import json
     import numpy as np
     from colorspace.colorlib import HCL, hexcols, colorobject
     from typing import List

     def generate_series(base:colorobject, n:int) -> List[str]:
         base.to('HCL')
         hues : np.ndarray = np.linspace(0,360, n+1)+base.get('H')[0]
         hues[hues >= 360] -= 360

         colors = HCL(H = hues[:-1], C = [base.get('C')[0]]*n, L = [base.get('L')[0]]*n)
         colors.to('hex')
         return colors.colors()

     base = hexcols('#AED882')
     colors = {i: generate_series(base, i) for i in range(1, 101)}
     with open('divergent_colors.json', 'w') as fh:
         json.dump(colors, fh)
"""

__all__ = ['feijoa', 'divergent_colors']

feijoa = '#AED882'

import json
from pathlib import Path
from typing import List, Dict

json_str = (Path(__file__).parent / 'divergent_colors.json').read_text(encoding='utf-8')

divergent_colors:Dict[int, List[str]] = {int(k): v for k, v in json.loads(json_str).items()}
