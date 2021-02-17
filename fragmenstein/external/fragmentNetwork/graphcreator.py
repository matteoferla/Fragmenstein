'''
Modified from https://github.com/xchem/fragalysis-api/blob/master/fragalysis_api/xcglobalscripts/config.ini and
https://github.com/xchem/fragalysis-api/blob/8946077b5cb20e4eb5409414fcc5ad61bbf222ac/fragalysis_api/xcanalyser/graphcreator.py and
https://github.com/xchem/strucbio_practical/blob/master/interactive_exercises/2-Compound_Design_Elaboration.ipynb
'''
import urllib.request as urllib
import json
from urllib.error import HTTPError

import pandas as pd

try:
    from fragalysis_api import ConfigSetup
    settings = ConfigSetup()
except ImportError:
    class Settings():
        _settings = {
            'fragalysis':{"url":"https://fragalysis.diamond.ac.uk/"},
            'graph':{"search":"network/full_graph/",
                     "query":"?smiles="}
        }
        @classmethod
        def get(cls, *args):
            current = Settings._settings
            for arg in args:
                current = current[arg]
            return current

    settings = Settings()


def xcgraphcreator(target_smiles):
    search = GraphRequest()

    new_smiles = search.get_new_smiles(smiles=target_smiles)

    return new_smiles


class GraphRequest:
    def __init__(self):

        # get url pieces
        self.frag_url = settings.get('fragalysis', 'url')
        self.graph_url = settings.get('graph', 'search')
        self.query = settings.get('graph', 'query')

        # get full url
        self.search_url = str(self.frag_url + self.graph_url + self.query)

        # set blanks for smiles search and json to handle later
        self.smiles_url = None
        self.graph_json = None

    def set_smiles_url(self, smiles):
        # set full search url
        self.smiles_url = str(self.search_url + smiles)

    def get_new_smiles(self, smiles):
        # check for a smiles url
        smiles_url = str(self.search_url + smiles)
        if not smiles_url:
            raise Exception('Please initiate smiles url with set_smiles_url(<smiles>)!')

        # get response from url and decode -> json
        with urllib.urlopen(smiles_url) as f:
            result = f.read().decode('utf-8')
            if not result == 'EMPTY RESULT SET':
                response = json.loads(result)

                # set json as decoded response for processing
                return graph_dict_to_df(response).reset_index(drop=True).copy()


# to flatten into a list for processing
def flatten_json(y):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)

    return out


def graph_dict_to_df(graph_dict):
    """
    This is the staircase to heaven

    :param graph_dict:
    :return:
    """
    a_df = pd.DataFrame()
    columns = ['type', 'insert_smiles', 'new_smiles', 'insertion']

    start = '2'

    for i1 in graph_dict[start].keys():
        for i2 in graph_dict[start][i1].keys():
            if isinstance(graph_dict[start][i1][i2], dict):
                for i3 in graph_dict[start][i1][i2].keys():
                    if isinstance(graph_dict[start][i1][i2][i3], dict):
                        for i4 in graph_dict[start][i1][i2][i3].keys():
                            if isinstance(graph_dict[start][i1][i2][i3][i4], list):
                                for i5 in graph_dict[start][i1][i2][i3][i4]:
                                    tmp_df = pd.DataFrame([[i1, i2, i5['end'], i5['change']]], columns=columns)
                                    a_df = pd.concat([a_df, tmp_df])

    return a_df


if __name__ == "__main__":
    import sys, os
    smiles = os.path.expanduser( sys.argv[1] )
    graph_search = xcgraphcreator(smiles)
    print(graph_search)
    print( set(graph_search.type.values) )

    '''
python -m fragmenstein.external.fragmentNetwork.graphcreator "CC1=NCCC(C)O1"

    '''