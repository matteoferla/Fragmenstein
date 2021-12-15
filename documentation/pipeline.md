> This code has not been updated. So usage may differ.

Here is an example of a pipeline to iterate across all the compounds in a list, merge them, score then and upload them.
For MPro, see the MPro.


## Prerequisite

### Template
I made an energy minimised template.
Technically, it need not be substrate bound, but having a closed active site is good.
Also it actually does not need to be energy minimised.

    pdbcode = '6WOJ' # ADP bound.
    
    from fragmenstein import Victor, Igor
    Igor.download_map(pdbcode, pdbcode+'.ccp4')
    
    # topology
    from rdkit_to_params import Params

    p = Params.from_smiles_w_pdbfile(pdb_file='mono.pdb', 
                                 smiles='Nc1ncnc2n(cnc12)[C@@H]3O[C@H](CO[P@@]([O-])(=O)O[P@@]([O-])(=O)OC[C@H]4O[C@@H](O)[C@H](O)[C@@H]4O)[C@@H](O)[C@H]3O',
                                 name='APR',
                                 proximityBonding=False)
    p.dump('APR.params')
    
    # start
    import pyrosetta
    pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
    
    #import nglview
    #nglview.show_rosetta(p.test())
    
    params_file = 'APR.params'
    pdbfile = 'mono.pdb' # 6WOJ manually inspected.
    
    pose = pyrosetta.Pose()
    params_paths = pyrosetta.rosetta.utility.vector1_string()
    params_paths.extend([params_file])
    pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdbfile)
    Igor.relax_with_ED(pose=pose, ccp4_file=pdbcode+'.ccp4')
    pose.dump_pdb('mono.r.pdb')
    
    # this part makes no sense, but is just an example —in reality checking the result visually would make more sense.
    # or it could be done in Pyrosetta.
    
    import pymol2
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load('mono.r.pdb')
        pymol.cmd.remove('resn APR')
        pymol.cmd.save('template.pdb')
    
### Compound extraction
XChem did not used to provide extracted compounds. These compounds below however preserve covalent bonding
by using the */R atom as the attachment. Say `CC[Cl]` reacted with `[S-]CC`, the former would be `CC*`, where `*` has the place of `S`.

    import os
    from rdkit import Chem
    from shutil import copyfile
    
    hits = {}
    # os.mkdir('comprimenda')
    masterfolder = '/Users/matteo/Desktop/NSP3-macrodomain/mArh/aligned'
    for subfolder in os.listdir(masterfolder):
        molfile = os.path.join(masterfolder, subfolder, f'{subfolder}.mol')
        if os.path.exists(molfile):
            mol = Chem.MolFromMolFile(molfile)
            mol.SetProp('_Name', subfolder)
            hits[subfolder] = mol
            copyfile(molfile, os.path.join('comprimenda', f'{subfolder}.mol'))


## Compound merger generation
It runs on a node on N processes and save the data as it goes along to a sqlite file to prevent and issues.

    project = 'mergers'
    nput_foldername = 'input'
    ##############################################
    cores = 20
    out_path = f'{project}'
    db_name = f'{project}.sqlite'
    ##############################################

Touch the directory and DB file

    import os, re
    from sqlitedict import SqliteDict
    import json
    results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)

definite the process task

    def process(x):
        project = 'mergers'
        db_name = f'{project}.sqlite'
        import pyrosetta, logging
        pyrosetta.distributed.maybe_init(extra_options='-no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false')
        from fragmenstein import Victor
        Victor.work_path = f'{project}'  # db_name
        Victor.fragmenstein_throw_on_discard= True
        Victor.fragmenstein_joining_cutoff = 5 # 10
        Victor.quick_renanimation = False
        Victor.error_to_catch = Exception
        #Victor.enable_stdout(logging.ERROR)
        Victor.enable_logfile(f'{project}.log', logging.INFO)
        Victor.log_errors()
        from sqlitedict import SqliteDict
        import json, logging
        from fragmenstein.mpro import MProVictor
        print('NEW', x)
        try:
            from rdkit import Chem
            
            def loadmol(file):
                mol = Chem.MolFromMolFile(file)
                if mol.GetProp('_Name') == '':
                    mol.SetProp('_Name', file.split('/')[-1].replace('.mol',''))
                return mol
            
            frags = [loadmol(file) for file in x]
            v = Victor(pdb_filename='input/template.pdb',
                        covalent_resi='81A', # a random residue is still required for the constaint ref atom.
                        covalent_resn='CYS')
            v.combine(hits=frags)
            results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
            results[v.long_name] = v.summarise()
            if not v.error:
                v.make_pse()
        except Exception as error:
            name = '-'.join([file.split('/')[-1].replace('.mol','') for file in x])
            error_msg = f'{error.__class__.__name__} {error}'
            results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
            name = '-'.join([file.split('/')[-1].replace('.mol','') for file in x])
            results[name] = {'error': error_msg}
            Victor.journal.critical(f'*** {error_msg}, files: {x}')
        except ConnectionError:
            pass
        print('DONE', x)
        return True
        
Get stuff started

    from multiprocessing import Pool
    import itertools, random, re
    pool = Pool(cores, maxtasksperchild=1)
    
    # get done list to prevent repeats
    results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
    done = list(results.keys())

    hits = [os.path.join(input_foldername, file) for file in os.listdir(input_foldername) if '.mol' in file]
    to_do = [(a, b) for a, b in itertools.permutations(hits, 2) if f'{a.split("/")[-1]}-{b.split("/")[-1]}'.replace('.mol', '') not in done]
    random.shuffle(to_do)
    for pair in to_do:
        pool.apply_async(process, (pair,))
        
The main process will be free. But the pool will be running.

    pool.close()
    pool.join()
    print("completed")
        
## Tabulate
The upload step requires, `michelanglo_api` uploading the data to github.
Method to make table:

    from sqlitedict import SqliteDict
    from rdkit.Chem import PandasTools
    import json
    import pandas as pd
    from fragmenstein import Victor
    
    import numpy as np
    import os
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from scipy.stats import skewnorm, gennorm
    
    
    def old_ranker(row):
        try:
            return float(row['∆∆G'])/5 + float(row.comRMSD) + row.N_unconstrained_atoms /5 - row.N_constrained_atoms / 10
            #return float(row['∆∆G'])/(row.N_unconstrained_atoms + row.N_constrained_atoms * 0.5)*10 + float(row.comRMSD)
        except:
            return float('nan')
        
    
    rank_weights = {'LE': 1., 'comRMSD': 2., 'atom_bonus': 2. , 'novelty_penalty': 5.}
    def ranker(row):
        try:
            #atom_bonus = row.N_constrained_atoms / (20 + row.N_constrained_atoms)
            #atom_bonus = skewnorm.pdf((row.N_constrained_atoms - 20)/8, 3)
            ζ = (row.N_constrained_atoms**2 - 25**2)/500
            atom_bonus = gennorm.pdf(ζ, 5) / 0.5445622105291682
            novelty_penalty = row.N_unconstrained_atoms / row.N_constrained_atoms
            return rank_weights['LE'] * float(row.LE) + \
                   rank_weights['comRMSD'] * float(row.comRMSD) + \
                   - rank_weights['atom_bonus'] * atom_bonus + \
                    rank_weights['novelty_penalty'] * novelty_penalty
        except:
            return float('nan')
        
    def LE(row):
        try:
            return float(row['∆∆G'])/(row.N_unconstrained_atoms + row.N_constrained_atoms)
        except:
            return float('nan')
    
    def get_mol3D(name):
        path = os.path.join(Victor.work_path, name, name+'.minimised.mol')
        if os.path.exists(path):
            try:
                mol = Chem.MolFromMolFile(path, sanitize=True)
                if mol is None:
                    return None
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                return mol
            except Exception as error:
                print(f'{type(error)}: {error}')
                return None
        else:
            return None
    
    
    def get_table(db_name, mols=True, mol_only=True):
        results = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
        result_table = pd.DataFrame(results.values())
        print(len(result_table), sum(~result_table['∆∆G'].isna()))
        result_table['LE'] = result_table.apply(LE,1)
        rank = result_table.apply(ranker, axis=1).rank()
        m = np.nanmax(rank.values)
        result_table['%Rank'] = rank / m * 100
        result_table['N_hits'] = result_table.regarded.apply(lambda x: len(x) if str(x) != 'nan' else float('nan'))
        result_table = result_table.loc[~result_table.smiles.isna()].sort_values(['%Rank'], axis=0) 
        if mols:
            result_table['mol3D'] = result_table['name'].apply(get_mol3D)
            #result_table['mol2D'] = result_table['name'].apply(get_mol2D)
            PandasTools.AddMoleculeColumnToFrame(result_table,'smiles','mol2D')
            if mol_only:
                result_table = result_table.loc[~result_table.mol3D.isna()]
        return result_table

Make table

    ##############################################
    project = 'mergers'
    from fragmenstein import Victor
    Victor.work_path = project
    db_name = f'{project}.sqlite'
    result_table = get_table(db_name, mols=True)
    #result_table
    
    # shoot. I forgot to count atoms in hits...
    atom_Ns = {}
    for folder in ('newinputs',): #'input', 'UCSF2-hits', 'frags'):
        for file in os.listdir(folder):
            if '.mol' in file:
                mol = Chem.MolFromMolFile(os.path.join(folder, file), sanitize=False)
                if mol is None:
                    atom_Ns[file.replace('.mol','')] = 0 # float nan?
                else:
                    mol = Chem.GetMolFrags(mol, asMols=True)[0] # just in case
                    atom_Ns[file.replace('.mol','')] = mol.GetNumAtoms()
                    
    # add atom n
    hit_counter = lambda hits: sum([atom_Ns[hit] for hit in hits])
    merge_counter = lambda row: row.N_hit_atoms - row.N_unconstrained_atoms - row.N_constrained_atoms 
    result_table = result_table.assign(N_hit_atoms=result_table.regarded.apply(hit_counter))
    result_table = result_table.assign(N_diff_atoms=result_table.apply(merge_counter, axis='columns'))
    result_table
    
    # upload
    from michelanglo_api import MikeAPI
    mike = MikeAPI('xxx', 'xxxx')
    # p = mike.convert_pdb('6WOJ') # make new
    p = mike.get_page('xxxxxxxx')  # retrieve old
    p.retrieve()
    p.show_link()

    task = 'Fragmenstein_NSP3_mergers'
    #repo_name = 'Data_for_own_Michelanglo_pages2'
    repo_name = 'NSP3-macrodomain'
    folder = 'molfiles'
    title = 'Fragmenstein NSP3 Mergers'
    #gitfolder=f'/well/brc/matteo/{repo_name}'
    gitfolder=f'/well/brc/matteo/NSP3/git-repo'
    sdfile=f'{gitfolder}/{folder}/mergers.sdf'
    xlsfile=f'{gitfolder}/{folder}/mergers.xlsx'
    targetfolder=f'{gitfolder}/{folder}'
    
    import os, re
    if not os.path.exists(targetfolder):
        os.mkdir(targetfolder)
        
    target = 'Mike'
#target = 'Excel'
#target = 'Frag'

    def clean_refs(x):
        h = ','.join([xx for xx in x if 'diamond-x' in xx])
        #h = ','.join(x)
        if h == '':
            return 'x0104_0'
        else:
            return h
    
    from datetime import datetime, date
    
    outgoing = result_table.sort_values(['%Rank'], axis=0).loc[~result_table.mol3D.isna()]
    outgoing = outgoing.loc[outgoing.regarded.apply(lambda x: len(x) >= 1)]
    outgoing['ref_mols'] = outgoing.regarded.apply(clean_refs)
    # frgaments have _Greek
    outgoing.ref_mols.replace('_[α-ω]', '', regex=True, inplace=True)
    outgoing.ref_mols.replace('mArh-', '', regex=True, inplace=True)
    outgoing.ref_mols.replace('diamond-', '', regex=True, inplace=True)
    outgoing.ref_mols.replace('mac-x(\d+)', r'MAC-\1_0_A', regex=True, inplace=True)
    outgoing['regarded'] = outgoing.regarded.apply(lambda x: ','.join(x))
    outgoing['disregarded'] = outgoing.disregarded.apply(lambda x: ','.join(x))
    # outgoing['ref_mols'] = outgoing.apply(lambda row: row.regarded+','+row.disregarded if row.disregarded else row.regarded, axis=1)
    outgoing['name'] = outgoing['name'].str.replace('-covalent', '')
    outgoing['name'] = outgoing['name'].str.replace('mArh-', '')
    outgoing['ref_pdb'] = 'x0104_0' #'6WOJ_0_A' #'X0104_0_A' #
    outgoing['original SMILES'] = outgoing['smiles']
    # ref_url - the url to the forum post that describes the work
    # submitter_name - the name of the person submitting the compounds
    # submitter_email - the email address of the submitter
    # submitter_institution - the submitters institution
    # generation_date - the date that the file was generated in format yyyy-mm-dd
    # method
    #ref_mols - a comma separated list of the fragments that inspired the design of the new molecule (codes as they appear in fragalysis - e.g. x0104_0,x0692_0)
    #ref_pdb - either (a) a filepath (relative to the sdf file) to an uploaded pdb file (e.g. Mpro-x0692_0/Mpro-x0692_0_apo.pdb) or (b) the code to the fragment pdb from fragalysis that should be used (e.g. x0692_0)
    #original SMILES
    
    
    metadata = {'name':           'ver_1.2',
                'submitter_name': 'Matteo Ferla',
                'submitter_email': 'matteo@well.ox.ac.uk',
                'submitter_institution': 'Universtity of Oxford',
                'generation_date': date.today().isoformat(),
                'method': task, #'Fragmenstein-top500-automerger',
                'original SMILES': 'molecule smiles',
                #'smiles': 'molecule smiles used',
                'ref_url': 'https://github.com/matteoferla/NSP3-macrodomain',
                'ref_mols': 'all reference molecules',
                'ref_pdb': 'All ligands were evaluated against the apo of Rosetta-ED-guided-minimised 6WOJ',
                'N_hits': 'Number of hits used',
                'N_constrained_atoms': 'Number of atoms in the submission that were constrained',
                'N_diff_atoms': 'Difference in number of heavy atoms between the merger and the hits (negative: atoms added, positive: atoms merged)',
                #'N_unconstrained_atoms': 'Number of heavy atoms in the submission that were NOT constrained',
                #'runtime': 'seconds it took',
                'regarded': 'Fragments used for mapping',
                'disregarded': 'Fragments rejected for mapping',
                'comRMSD': 'Combined RMSD from the atoms of the fragments that contributed to the position of the followup',
                '∆∆G': 'Difference in Gibbs Free energy relative to unbound molecule in kcal/mol (ref2015 scorefxn; negative=Good)',
                #'∆G_bound': 'Gibbs Free energy of ligand bound',
                #'∆G_unbound': 'Gibbs Free energy of ligand unbound',
                'LE': 'Ligand efficiency (kcal/mol/N_heavy)',
                '%Rank':f"Sorted by {rank_weights['comRMSD']}x RSMD (high is bad) + "+\
                        f"{rank_weights['LE']}x ligand efficiency (high is bad) - "+\
                        f"{rank_weights['atom_bonus']}x N_constrained_atoms/(20+N_constrained_atoms) + "+\
                        f"{rank_weights['novelty_penalty']}x N_unconstrained_atoms/N_constrained_atoms",
                'mol3D': Chem.MolFromSmiles('FOO'),
                'mol2D': Chem.MolFromSmiles('FOO'),
              }
    
    if target == 'Frag':
        del metadata['regarded']
        del metadata['disregarded']
        del metadata['mol2D']
        outgoing = outgoing.iloc[:500]
    elif target == 'Mike':
        del metadata['mol2D']
    elif target == 'Excel':
        del metadata['ref_url']
        del metadata['ref_mols']
        del metadata['ref_pdb']
        del metadata['submitter_name']
        del metadata['submitter_email']
        del metadata['submitter_institution']
        del metadata['generation_date']
        del metadata['mol3D']
        
    #outgoing = outgoing[list(set(outgoing.columns.values) - set(metadata.keys()))]
    # add fist compound
    outgoing = pd.concat([pd.DataFrame([metadata]), outgoing], ignore_index=True)

    # upload
    title = 'Fragmenstein/Smallworld NSP3 mergers'

    p.description = f'''
    ## {title} {date.today().isoformat()}
    
    [Fragmenstein](https://github.com/matteoferla/Fragmenstein) scored suggestions from Smallworld server.
    
    > For files and notes see [NSP3 data on GitHub](https://github.com/matteoferla/NSP3-macrodomain).
    
    '''
    p.loadfun = ''
    p.title = title
    p.columns_viewport = 6
    p.columns_text = 6
    p.sdf_to_mols(sdfile=sdfile,
                     targetfolder=targetfolder,
                     skip_first=True)
    
    
    
    p.sdf_to_json(sdfile=sdfile,
                     keys=('regarded',
                           '∆∆G', 'LE', 'N_hits', 'N_constrained_atoms', 'N_diff_atoms', 'comRMSD',
                           '%Rank'),
                     key_defaults=('', # regarded
                                   999., #∆∆G
                                   999., #LE
                                   0, #N_hits
                                   0, #N_constrained_atoms
                                   0, #N_diff_atoms
                                   999., #comRMSD
                                   100. #%Rank
                                  ), 
                     filename=f'{targetfolder}/data.json',
                     spaced=True)
    p.make_fragment_table(sdfile=sdfile,
                       username='matteoferla',
                       repo_name=repo_name,
                       foldername=folder,
                       protein_sele='81:A',
                       sort_col=8,
                       sort_dir='asc',
                       template_row=-1,
                       fragment_row=1,
                       jsonfile='data.json')
    p.commit()

Copy the inputs

    import os, re, shutil
    
    hit_folder = 'newinputs' #'UCSF2-hits' #'frags' #'input'
    
    for file in os.listdir(hit_folder):
        if '.mol' in file:
            print(file)
            shutil.copy(os.path.join(hit_folder, file),
                       os.path.join(targetfolder, file.replace('_0_','_'))) #
    shutil.copy(os.path.join('input', 'template.pdb'),
                       os.path.join(targetfolder, 'template.pdb'))

Git push!

Alternatively for an SDF for Fragalysis say

    ## MAKE SDF
    assert target != 'Excel', 'Requires mol2D'
    from rdkit.Chem import PandasTools
    #Fragmenstein_permissive_rescored_20200609.sdf
    PandasTools.WriteSDF(outgoing.iloc[:501], sdfile, molColName='mol3D', idName='name',
                         properties=list(set(metadata.keys()) - {'name', 'mol3D', 
    #                                                              'N_diff_atoms', 
                                                                 'method', 
                                                                 'submitter_name',
                                                                'submitter_institution',
                                                                'submitter_email'
                                                                }), allNumeric=False)