from typing import Dict, List, Union


def process(data: Dict[str, Union[str, dict]]):
    # read data
    project = data['project']
    db_name = f'{project}.sqlite'
    hit_blocks = data['hit_blocks']
    name = '-'.join(hit_blocks.keys())
    print('**********', name)
    # imports ------------------------------------------------------
    import pyrosetta, logging
    from .make_pyrosetta_options import make_option_string
    pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                                      mute='all',
                                                                      ignore_unrecognized_res=True,
                                                                      load_PDB_components=False))
    from fragmenstein.mpro import MProVictor
    from sqlitedict import SqliteDict
    import json, logging
    # fix hits ------------------------------------------------------
    from rdkit import Chem
    hits = []
    for hit_name in hit_blocks:
        hit = Chem.MolFromMolBlock(hit_blocks[hit_name])
        hit.SetProp('_Name', hit_name)
        hits.append(hit)
    # settings for Fragmenstein ------------------------------------
    MProVictor.work_path = f'{project}'  # db_name
    MProVictor.monster_throw_on_discard = True
    MProVictor.joining_cutoff = 5  # 10
    MProVictor.quick_reanimation = False
    MProVictor.error_to_catch = Exception
    MProVictor.enable_stdout(logging.ERROR)
    MProVictor.enable_logfile(f'{project}.log', logging.INFO)
    MProVictor.capture_rdkit_log()
    # analyse ------------------------------------------------------
    try:
        v = MProVictor.combine(hits=hits)
        results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
        results[v.long_name] = v.summarize()
        if not v.error_msg:
            v.make_pse()
        print('DONE', [hit.GetProp('_Name') for hit in hits])
        return v.minimised_mol
    except Exception as error:
        error_msg = f'{error.__class__.__name__} {error}'
        results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
        results[name] = {'error': error_msg}
        MProVictor.journal.critical(f'*** {error_msg}, data: {data}')
    except ConnectionError:
        pass
    return None
