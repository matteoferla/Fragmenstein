import json, gzip
import numpy as np

from fragmenstein.external import ExternalToolImporter

[standalone_model_numpy] = ExternalToolImporter.import_tool("scscore", ["scscore.standalone_model_numpy"])


class SCScoreWrapper():
    '''
    A wrapper class to use SCScore in a transparent manner.
    '''



    def __init__(self, ):
        '''
        '''

        def _load_vars(other, weight_path):
            if weight_path.endswith('pickle'):
                import pickle
                with open(weight_path, 'rb') as fid:
                    other.vars = pickle.load(fid, encoding= 'latin1')
                    other.vars = [x.tolist() for x in other.vars]
            elif weight_path.endswith('json.gz'):
                with gzip.GzipFile(weight_path, 'r') as fin:  # 4. gzip
                    json_bytes = fin.read()  # 3. bytes (i.e. UTF-8)
                    json_str = json_bytes.decode('utf-8')  # 2. string (i.e. JSON)
                    other.vars = json.loads(json_str)
                    other.vars = [np.array(x) for x in other.vars]

        standalone_model_numpy.SCScorer._load_vars = _load_vars

        self.scscorer = standalone_model_numpy.SCScorer()

        self.scscorer.restore()

    def compute_score(self, *smiles):
        scores = []
        for smi in smiles:
            score = self.scscorer.get_score_from_smi(smi)[-1]
            scores.append( score )

        return scores



if __name__ == "__main__":
    import sys
    scsw = SCScoreWrapper()
    scores = scsw.compute_score( sys.argv[1:] )
    print(scores)

    '''
python -m fragmenstein.external.scscore.SCScoreWrapper
    '''