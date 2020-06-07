import pyrosetta

pyrosetta.init(extra_options='-no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false')

from fragmenstein.mpro import MProVictor

MProVictor.work_path = 'mpro_output'

from sqlitedict import SqliteDict

import json

import logging

MProVictor.enable_stdout(logging.DEBUG)


#####################################################
results = SqliteDict('results.sqlite', encode=json.dumps, decode=json.loads, autocommit=True)

# data from https://github.com/postera-ai/COVID_moonshot_submissions
postera = MProVictor.fetch_postera()

for i, row in postera.iterrows():
    if row.fragments == 'x0072':
        # these are not hit inspired.
        continue
    if row.covalent_warhead in (False, 'False', 'false'):
        v = MProVictor.from_hit_codes(long_name=row.CID,
                                      hit_codes=row.fragments.split(','),
                                      smiles=row.SMILES)
        v.make_pse()
        results[row.CID] = v.summarise()
    else:
        if row.category not in ('Acrylamide', 'Chloroacetamide', 'Vinylsulfonamide', 'Nitrile'):
            print(f'What is {row["CID"]}')
        else:
            combinations = MProVictor.make_all_warhead_combinations(row.SMILES, row.category)
            if combinations is None:
                break
            for c in combinations:
                if '_noncovalent' in c:
                    pass  # the pre-encounter complex needs to be done differently (Code to be done)
                else:
                    v = MProVictor.from_hit_codes(long_name=row.CID + '-' + c,
                                                  hit_codes=row.fragments.split(','),
                                                  smiles=combinations[c])
                    v.make_pse()
                    results[row.CID] = v.summarise()

# requires slack webhook. See api.slack.com
MProVictor.slack_me('DONE')