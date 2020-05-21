import pyrosetta

pyrosetta.init(extra_options='-no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false')  #
from fragmenstein import Igor, Fragmenstein, Victor
import logging, csv, json, random
from rdkit import Chem
from rdkit.Chem import AllChem
from sqlitedict import SqliteDict

results = SqliteDict('results.sqlite', encode=json.dumps, decode=json.loads, autocommit=True)
########################################################################################################################

Victor.fragmenstein_merging_mode = 'none'
Victor.work_path = 'mpro_output'
# Victor.enable_stdout(logging.DEBUG)
Victor.enable_logfile('reanimate.log')

for cname, con in [('chloroacetamide', 'AtomPair  H  145A  OY  1B HARMONIC 2.1 0.2\n'),
                   ('nitrile', 'AtomPair  H  145A  NX  1B HARMONIC 2.1 0.2\n'),
                   ('acrylamide', 'AtomPair  H  143A  OZ  1B HARMONIC 2.1 0.2\n'),
                   ('vinylsulfonamide', 'AtomPair  H  143A  OZ1 1B HARMONIC 2.1 0.2\n')
                   ]:
    Victor.add_constraint_to_warhead(name=cname, constraint=con)

mpro_folder = 'test_data/Mpro_mols_18May'


def get_category(row):
    for category in ('acrylamide', 'chloroacetamide', 'vinylsulfonamide', 'nitrile'):
        if row[category] == 'True':
            return category
    else:
        return None


def get_mol(xnumber):
    xnumber = xnumber.strip()
    mol = Chem.MolFromMolFile(f'{mpro_folder}/Mpro-{xnumber}_0/Mpro-{xnumber}_0.mol')
    mol.SetProp('_Name', xnumber)
    return mol


# def get_best(hit_codes):
#     return Victor.closest_hit(pdb_filenames=[f'{mpro_folder}/Mpro-{i}_0/Mpro-{i}_0_bound.pdb' for i in hit_codes],
#                               target_resi=145,
#                               target_chain='A',
#                               target_atomname='SG',
#                               ligand_resn='LIG')


def pose_fx(pose):
    pose2pdb = pose.pdb_info().pdb2pose
    r = pose2pdb(res=41, chain='A')
    MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
    MutateResidue(target=r, new_res='HIS').apply(pose)


def poised_pose_fx(pose):
    pose2pdb = pose.pdb_info().pdb2pose
    r = pose2pdb(res=41, chain='A')
    MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
    MutateResidue(target=r, new_res='HIS_D').apply(pose)
    r = pose2pdb(res=145, chain='A')
    MutateResidue(target=r, new_res='CYZ').apply(pose)


def reanimate(smiles, name, hit_codes, category=None):
    hit_codes = [h for h in hit_codes if h not in ('x0748', 'x0970', 'x0947', 'x1420')]
    hits = [get_mol(i) for i in hit_codes]
    # best_hit = get_best(hit_codes)
    # Victor.journal.debug(f'{name} - best hit as starting is {best_hit}')
    # apo = best_hit.replace('_bound', '_apo-desolv')
    apo = 'template.pdb'
    atomnames = {}
    fx = pose_fx
    extra_constraint = 'AtomPair  SG  145A  NE2  41A HARMONIC 3.5 0.2\n'
    if category:
        cname, rxd = category.split('_')
        if rxd == 'noncovalent':
            wd = [wd for wd in Victor.warhead_definitions if wd['name'] == cname][0]
            mol = Chem.MolFromSmiles(smiles)
            nc = Chem.MolFromSmiles(wd['noncovalent'])
            atomnames = dict(zip(mol.GetSubstructMatch(nc), wd['noncovalent_atomnames']))
            fx = poised_pose_fx
            extra_constraint += 'AtomPair  SG  145A  CX   1B HARMONIC 3.2 0.5\n'
            extra_constraint += wd['constraint']
    print(f'reanimate(smiles="{smiles}", name="{name}", hit_codes={hit_codes})')
    reanimator = Victor(smiles=smiles,
                        hits=hits,
                        pdb_filename=apo,
                        long_name=name,
                        ligand_resn='LIG',
                        ligand_resi='1B',
                        covalent_resn='CYS', covalent_resi='145A',
                        extra_constraint=extra_constraint,
                        pose_fx=fx,
                        atomnames=atomnames
                        )
    return reanimator


def store(v):
    # persistent
    results[v.long_name] = {'name': v.long_name,
                            'mode': 'none',
                            '∆∆G': v.energy_score['ligand_ref2015']['total_score'] - \
                                   v.energy_score['unbound_ref2015']['total_score'],
                            'comRMSD': v.mrmsd.mrmsd,
                            'N_constrained_atoms': v.constrained_atoms,
                            'runtime': v.tock - v.tick,
                            'disregarded': json.loads(v.fragmenstein.scaffold.GetProp('parts'))
                            }


#####################################################

# data from https://github.com/postera-ai/COVID_moonshot_submissions
data = list(csv.DictReader(open('../COVID_moonshot_submissions/covid_submissions_all_info.csv')))

random.shuffle(data)

#######################################################

truths = (True, 'TRUE', 'True', 'true', 1, 'VERO', 'YES', 'Yes', 'yes', 'SI', 'Si', 'si', 'JA', 'Ja', 'Aye')
falsehoods = (False, 'FALSE', 'False', 'false', 'FALSO', 0, None, 'NO', 'No', 'no', 'NEIN', 'Nein', 'Nay')

for row in data:
    if row['fragments'] == 'x0072':
        continue
    if row['covalent_warhead'] in falsehoods:
        reanimate(name=row['CID'], hit_codes=row['fragments'].split(','), smiles=row['SMILES'])
    else:
        print(f'Covalent: {row["CID"]}')
        for category in ('acrylamide', 'chloroacetamide', 'vinylsulfonamide', 'nitrile'):
            if row[category] == 'True':
                combinations = Victor.make_all_warhead_combinations(row['SMILES'], category)
                if combinations is None:
                    break
                for c in combinations:
                    if '_noncovalent' in c:
                        pass  # this needs to be done differently (Code to be done)
                        # reanimate(name = row['CID']+'-'+c+'_UNREACTED',
                        #           hit_codes = row['fragments'].split(','),
                        #           smiles=combinations[c],
                        #           category=c)
                    else:
                        reanimate(name=row['CID'] + '-' + c, hit_codes=row['fragments'].split(','),
                                  smiles=combinations[c])
                break
        else:
            print(f'What is {row["CID"]}')

# requires slack webhook. See api.slack.com
Victor.slack_me('DONE')
