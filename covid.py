##############################################
cores = 2
out_path = 'mpro_output'
##############################################

import pandas as pd
from multiprocessing import Pool


def preprocess(path):
    import pyrosetta, logging
    pyrosetta.init(extra_options='-no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false')
    from fragmenstein.mpro import MProVictor
    MProVictor.enable_stdout(logging.INFO)
    MProVictor.work_path = path

def process(x):
    from sqlitedict import SqliteDict
    import json
    from fragmenstein.mpro import MProVictor
    v = MProVictor.from_hit_codes(**x)
    v.make_pse()
    results = SqliteDict(f'MPro.sqlite', encode=json.dumps, decode=json.loads, autocommit=True)
    results[x['long_name']] = v.summarise()

def get_compounds_to_be_done(postera: pd.DataFrame):
    # parse
    to_be_done = []
    mystery = []
    for i, row in postera.iterrows():
        if row.fragments == 'x0072' or str(row.fragments) == 'nan':
            # these are not hit inspired.
            continue
        if row.covalent_warhead in (False, 'False', 'false'):
            # parse
            to_be_done.append(dict(long_name=row.CID,
                                          hit_codes=row.fragments.split(','),
                                          smiles=row.SMILES))
        else:
            if row.category.lower() not in ('acrylamide', 'chloroacetamide', 'vinylsulfonamide', 'nitrile'):
                print(f'What is {row["CID"]}')
                mystery.append(row["CID"])
            else:
                combinations = MProVictor.make_all_warhead_combinations(row.SMILES, row.category)
                if combinations is None:
                    mystery.append(row["CID"])
                    break
                for c in combinations:
                    if '_noncovalent' in c:
                        pass  # the pre-encounter complex needs to be done differently (Code to be done)
                    if row.category.lower() in c:
                        to_be_done.append(dict(long_name=row.CID,
                                               hit_codes=row.fragments.split(','),
                                               smiles=combinations[c]))
                    else:

                        to_be_done.append(dict(long_name=row.CID + '-' + c.split('_')[0],
                                                      hit_codes=row.fragments.split(','),
                                                      smiles=combinations[c]))
    print(mystery)
    return to_be_done

if __name__ == '__main__':
    # get stuff started
    pool = Pool(cores)
    #pool.map(preprocess, [out_path] * cores)

    # get data
    # data from https://github.com/postera-ai/COVID_moonshot_submissions
    from fragmenstein.mpro import MProVictor
    postera = MProVictor.fetch_postera()
    to_be_done = get_compounds_to_be_done(postera)
    # pool.map(process, to_be_done)
    # requires slack webhook. See api.slack.com
    MProVictor.slack_me('DONE')