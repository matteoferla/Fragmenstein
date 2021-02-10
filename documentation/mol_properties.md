## Molecule properties

> This file is for discussion of standards

### SDFile
As of February 2021, Monster stores information in the Properties fields of the molecule.
These are not saved as MOL files, but are in SDF files.

    import tempfile
    def _MolToSDBlock(mol:Chem.Mol):
        with tempfile.NamedTemporaryFile(suffix='.sdf') as fh:
            w = Chem.SDWriter(fh.name)
            return w.GetText(mol)
    Chem.MolToSDBlock = _MolToSDBlock

    mol = Chem.MolFromSmiles('*CCC')
    mol.SetProp('_Name','Something')
    mol.SetProp('_Hide','Nothing')
    mol.SetProp('Show','Everything')
    print('Mol file has no properties')
    print(Chem.MolToMolBlock(mol))
    print('SD file has non-private properties')
    print(Chem.MolToSDBlock(mol))
    
### Binary
A pickled Chem.Mol does not require sanitation, so for 
A further issue is that pickled `Chem.Mol` objects lack properties,
hence the presence of `PropertyMol` workaround ([link](https://www.rdkit.org/docs/source/rdkit.Chem.PropertyMol.html)).
However, the properties assigned to a `Mol` before conversion to a `PropertyMol` are lost.
So a better solution is to use the `ToBinary` method of a mol instance.

    bstr = mol.ToBinary(propertyFlags=0b00010111)
    with open('out.b', 'b') as fh:
        fh.write(bstr)
        
The binary string can be passed to a Chem.Mol:

    Chem.Mol(bstr)

I could not find any documentation for the exact definitions of property flags argument, but empirical I have figure out:

    0b00000001 -> mol
    0b00000010 -> atom
    0b00000100 -> bond
    0b00010001 -> inc. private
    
I am guessing that it's not a full octet of bits, but simply 5 bits. `0b1000` is probably include calculated.
I'd say this is probably the best way to store them, albeit non-standard.

### Custom
Alternative, one could have a custom extractor to save as JSON.


    def get_properties(mol: Chem.Mol) -> Dict[str, Union[dict, List[dict]]]:
        data = dict(mol= mol.GetPropsAsDict(includePrivate=True, includeComputed=False),
                    atoms=[atom.GetPropsAsDict(includePrivate=True, includeComputed=False) for atom in mol.GetAtoms()],
                    bonds=[bond.GetPropsAsDict(includePrivate=True, includeComputed=False) for bond in mol.GetBonds()])
        return data

    def set_properties(mol: Chem.Mol, data: Dict[str, dict]):
        def _assign(obj, k, v):
            if isinstance(v, str):
                obj.SetProp(k, v)
            elif isinstance(v, int):
                obj.SetIntProp(k, v)
            elif isinstance(v, bool):
                obj.SetBoolProp(k, v)
            elif isinstance(v, float):
                obj.SetDoubleProp(k, v)
            elif v is None:
                pass
            else:
                obj.SetProp(k, str(v))
        for k,v in data['mol'].items():
            _assign(mol, k, v)
        for atom, atomdata in zip(mol.GetAtoms(), data['atoms']):
            for k,v in atomdata.items():
                _assign(atom, k, v)
        for bond, bonddata in zip(mol.GetBonds(), data['bonds']):
            for k, v in bonddata.items():
                _assign(bond, k, v)