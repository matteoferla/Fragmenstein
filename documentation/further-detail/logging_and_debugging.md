## Logging

Logging is handled via the `logging` module. The handler is not defined however.
This is true for `Victor`, `Fragmenstein`, `Igor`, `Rectifier` and the dependency `Params`.

In `Victor` the logger is a class attribute called `.journal`, in the other it's simply `.log`.

Victor has a shorthand for adding a file or console handler.

    import logging
    Victor.enable_stdout(level=logging.DEBUG)
    Victor.enable_logfile(level=logging.DEBUG)

To change the level quickly just run it again.

## Error

If an error occurs, it will be caught and logged and the user will just get a log error.
If instead one wants to have the error be left uncaught:

    Victor.error_to_catch = ()
    
The error `ConnectionError` gets used internally (for bond problems, not internet ones) 
and is caught regardless of this.

## Steps

During merging and rectification, odd things are bound to happen.
Given a `Victor` instance `victor`, this could be done to make a pymol session to see what happened
at each step:

    import pymol2
    from rdkit import Chem
    with pymol2.PyMOL() as pymol:
        for hit in victor.fragmenstein.hits:
            # this may differ from victor.hits —the latter are not collapsed.
            pymol.cmd.read_molstr(Chem.MolToMolBlock(hit, kekulize=False), hit.GetProp('_Name'))
        for i, mod in enumerate(victor.modifications):
            pymol.cmd.read_molstr(Chem.MolToMolBlock(mod, kekulize=False), f'step{i}')
        pymol.cmd.save('test.pse')

The `Monster` class can trip up, and this can happen at various steps.
`victor.modifications` contains a stack of the various step, such as the collapsed rings to the rectification.


