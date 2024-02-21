## 3D Viewing

Two common 3D viewers running on different tech are `nglview` and `py3Dmol`.
Both are wrappers for JavaScript libraries and can be used in Jupyter notebooks.
The problem is that nglview looks nicer but currently breaks colab due to ipywidgets version.

Monster, Victor and Walton have a ``show``.
If either `nglview` or `py3Dmol` are installed, they will use them (cf. `fragmenstein.display.DISPLAYMODE`).
To specify which give the argument `viewer_mode` to show (`nglview` and `py3Dmol`).

`Monster.show()` will therefore call `Monster.to_nglview()` or `Monster.to_py3Dmol()` depending on the mode.
Both call `_add_mols_to_viewer` internally as the viewers are modified
(by subclassing in NGLWidget and by monkeypatching methods in py3Dmol) to behave similarly with their `add_mol` method.

`Victor.show()` just uses `Monster.show()` under the hood.

The `add_mol` has the attribute `carbon_color`.
The colours are based on "branding" by default. But can be overridden via Mol prop `color`.
py3Dmol has been altered to accept hex colours (cf. https://blog.matteoferla.com/2024/01/custom-carbon-colours-in-py3dmol.html)

In Walton there is a problem that if a compound is altered, the viewer is not updated.

However, given that this is a minor detail and py3Dmol does not work in a lot of cases, this will not be inestigated further.