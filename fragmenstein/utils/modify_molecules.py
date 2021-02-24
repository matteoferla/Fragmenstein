from rdkit import Chem


def remove_dummy_atom(mol, dummy_symbol= '*'):
    return Chem.DeleteSubstructs(mol, Chem.MolFromSmiles(dummy_symbol))

def change_atom_type(mol, initial_symbol= '*', final_symbol=["Xe"]):
    result = Chem.ReplaceSubstructs(mol, Chem.MolFromSmiles(initial_symbol), Chem.MolFromSmiles(final_symbol), replaceAll=True)[0]
    # print( result)
    return result
    # raise  NotImplementedError()