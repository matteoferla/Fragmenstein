## Ring collapse

A "simple" solution is to replace all rings with a single atom that can be unpacked later.

``Ring`` class in ``core._collapse_ring`` does exactly that (inherited by ``Frankenstein``).

![collapse](images/atom_collapse.png)

But is totally buggy in the upacking step after merging, therefore it is not implemented in Victor.

* `.collapse_ring(mol)`
* `.expand_ring(mol)`

Two problems I am working on slowly —this is an afterhours project:

* ~Double rings result in extra atoms at the ring!~ —fixed
* Ring collapsed merging results in some odd bonds!

Both have issues with `expand_ring`.

The following works now:

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

Whereas this does not:

    mpro_folder = '/Users/matteo/Coding/Mpro'
    
    def get_mol(xnumber):
        xnumber = xnumber.strip()
        mol = Chem.MolFromMolFile(f'{mpro_folder}/Mpro-{xnumber}_0/Mpro-{xnumber}_0.mol')
        mol.SetProp('_Name', xnumber)
        return mol
        
    a = get_mol('x0107')
    þ.store_positions(a)
    a = þ.collapse_ring(a)
    b = get_mol('x1093')
    þ.store_positions(b)
    b = þ.collapse_ring(b)
    m = þ.merge_pair(a, b)
    m2 = þ.expand_ring(m)
    
    import nglview
    view = nglview.NGLWidget()
    view.add_component(a)
    view.add_component(b)
    view
    
The problem is the bonding is not right.
This pair is nasty because one ring becomes a double one.
It requires the valence fix which is fine (as seen in postera site demo).
However, the bonds need to be fixed going to the new ring.
The match collapsed ring with collapsed ring only is clearly unneeded.
Proximity bonding is a nice way to resolve this.

## Logging

> Depracation of `Fragmenstein.notebook` (txt dict for debug) in favour of `Victor.journal` (proper logging)

The correct logging is via Victor's journal. The logging log `Fragmenstein`.
However, Fragmenstein has `.logbook` which are more debug notes —dictionary. These should be integrated or removed.

## Unittest

`test.py` and the various jupyter notebook I have need to go into a proper testing suite.