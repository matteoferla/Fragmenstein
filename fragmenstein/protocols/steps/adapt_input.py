from collections import OrderedDict

from fragmenstein.protocols.dataModel.compound import Compound


class InputAdapter():

    @classmethod
    def adapt_dict_or_compoundsList(cls, input_data ):
        assert  len(input_data) > 0 ,"Error, novalid input data provided"
        if isinstance(input_data, (tuple, list)) and isinstance(input_data[0], Compound):
            return OrderedDict([ (x.molId, x) for x in input_data])
        elif isinstance(input_data, dict):
            return OrderedDict([ (key, Compound(mol, molId=key)) for key, mol in input_data.items() ])
        else:
            raise ValueError("Error, not recognized input type %s : %s"%(input_data, type(input_data)))