from fragmenstein.utils.compound import Compound


class InputAdapter():
    def adapt_dict_or_compoundsList(self, input_data ):
        if isinstance(input_data, (tuple, list)) and isinstance(input_data[0], Compound):
            return {x.molId: x for x in input_data}
        elif isinstance(input_data, dict):
            return { key: Compound(mol, molId=key) for key, mol in input_data.items() }
        else:
            raise ValueError("Error, not recognized input type %s : %s"%(input_data, type(input_data)))