## Ring collapse

A "simple" solution is to replace all rings with a single atom that can be unpacked later.

``Ring`` class in ``core._collapse_ring`` does exactly that (inherited by ``Frankenstein``).

![collapse](images/atom_collapse.png)

But is totally buggy in the upacking step after merging, therefore it is not implemented in Victor.

* `.collapse_ring(mol)`
* `.expand_ring(mol)`

Two problems I am working on slowly —this is an afterhours project:

* Double rings result in extra atoms at the ring!
* Ring collapsed merging results in some odd bonds!

Both have issues with `expand_ring`.

    from fragmenstein.core import Ring
    from fragmenstein import Fragmenstein
    from IPython.display import display
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    class Þragmenstein(Fragmenstein): # mocked class.
            def __init__(self):
                self._debug_draw = False
    draw = Þragmenstein().draw_nicely
    
    for smiles in ['c1ccccc1CCCCN','c1ccccc1CCCCNCc2ccccc2', 'CCc1ccc2ccccc2c1']:
    m = Chem.MolFromSmiles(smiles)
    m.SetProp('_Name', smiles)
    AllChem.Compute2DCoords(m)
    ring = Ring(_debug_draw = True)
    ring.store_positions(m)
    print('original')
    display(m)
    mod = ring.collapse_ring(m)
    print('collapsed')
    display(mod)
    print('expanded')
    mod2 = ring.expand_ring(mod)
    display(mod2)
