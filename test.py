# These are not unit-tests!
import unittest, os


import pyrosetta

pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

from fragmenstein import Monster, Victor, Igor, Rectifier
from fragmenstein.mpro import MProVictor
from rdkit import Chem
from rdkit.Chem import AllChem
import logging
Victor.enable_stdout(level=logging.DEBUG)
Victor.work_path = 'test_output'


def test_molecule(name, smiles, hitnames):
    hits = [Chem.MolFromMolFile(f'../Mpro/Mpro-{i}_0/Mpro-{i}_0.mol') for i in hitnames]
    followup = Chem.MolFromSmiles(smiles)
    r = Chem.MolFromMolFile(f'../Mpro/Mpro-{hitnames[0]}_0/Mpro-{hitnames[0]}_0_SG.mol')
    f = Monster(followup, hits, attachment=r)
    f.make_pse(f'test_{name}.pse')


def easy_test():
    test_molecule(name='2_ACL', smiles='CCNc1ncc(C#N)cc1CN1CCN(C(=O)C*)CC1', hitnames=('x0692', 'x0305', 'x1249'))


def nasty_test():
    test_molecule(name='AGN-NEW-5f0-1_ACR1',
                  smiles='*CCC(=O)N1CC(CCN(C(=O)Nc2c(C)ncc(C)c2CCN2CCOCC2)c2cc(C)ccn2)C1',
                  hitnames=('x0434', 'x0540'))


def nasty2_test():
    test_molecule(name='BEN-VAN-c98-4',
                  smiles='*C(=N)CN1CCN(Cc2ccc(-c3cc(CC)ncn3)c(F)c2)CC1',
                  hitnames='x0692,x0770,x0995'.split(','))

def igor_test():
        acl = Igor.from_pdbfile(pdbfile='output/PAU-WEI-b9b-8_NIT2/pre_PAU-WEI-b9b-8_NIT2.pdb',
                                params_file='output/PAU-WEI-b9b-8_NIT2/PAU-WEI-b9b-8_NIT2.params',
                                constraint_file='cysBound.noCY.cst')
        acl.coordinate_constraint = 100
        print('initial')
        print(acl.ligand_score())
        print(acl.ligand_residue)  # AtomPair SG 145A HE2 41A HARMONIC 1.8 0.2
        print(acl.key_residues)
        r = acl.pose.pdb_info().pdb2pose(res=41, chain='A')
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res='HIS_D').apply(acl.pose)
        MutateResidue(target=r, new_res='HIS').apply(acl.pose)
        pymover = pyrosetta.PyMOLMover()
        pymover.pymol_name(f'initial')
        native = acl.pose.clone()
        pymover.apply(native)
        acl.pose = native.clone()
        acl.minimise()
        pymover.pymol_name(f'relaxed')
        pymover.apply(acl.pose)
        print('cartesian relaxed')
        print(acl.ligand_score())
        acl.pose.dump_pdb('igor_test.pdb')

def victor_test():
    Victor.work_path = '../Mpro_fragmenstein'
    Victor.enable_stdout(logging.DEBUG)

    for cname, con in [('chloroacetamide', 'AtomPair H 145A OY 1B HARMONIC 2.1 0.2\n'),
                       ('nitrile', 'AtomPair H 145A NX 1B HARMONIC 2.1 0.2\n'),
                       ('acrylamide', 'AtomPair H 143A OZ 1B HARMONIC 2.1 0.2\n'),
                       ('vinylsulfonamide', 'AtomPair H 143A OZ1 1B HARMONIC 2.1 0.2\n')
                       ]:
        Victor.add_constraint_to_warhead(name=cname, constraint=con)

    mpro_folder = '/Users/matteo/Coding/rosettaOps/Mpro'

    def get_mol(xnumber):
        mol = Chem.MolFromMolFile(f'{mpro_folder}/Mpro-{xnumber}_0/Mpro-{xnumber}_0.mol')
        mol.SetProp('_Name', xnumber)
        return mol

    def get_best(hit_codes):
        return Victor.closest_hit(pdb_filenames=[f'{mpro_folder}/Mpro-{i}_0/Mpro-{i}_0_bound.pdb' for i in hit_codes],
                                  target_resi=145,
                                  target_chain='A',
                                  target_atomname='SG',
                                  ligand_resn='LIG')

    def pose_fx(pose):
        pose2pdb = pose.pdb_info().pdb2pose
        r = pose2pdb(res=41, chain='A')
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res='HIS').apply(pose)

    def reanimate(smiles, name, hit_codes):
        hits = [get_mol(i) for i in hit_codes]
        best_hit = get_best(hit_codes)
        Victor.journal.debug(f'{name} - best hit as starting is {best_hit}')
        apo = best_hit.replace('_bound', '_apo-desolv')
        reanimator = Victor(smiles=smiles,
                            hits=hits,
                            pdb_filename=apo,
                            long_name=name,
                            ligand_resn='LIG',
                            ligand_resi='1B',
                            covalent_resn='CYS', covalent_resi='145A',
                            extra_constraint='AtomPair SG 145A NE2 41A HARMONIC 3.5 0.2\n',
                            pose_fx=pose_fx
                            )
        return reanimator

    reanimate(name='DAV-CRI-d1e-2_ACR',
              hit_codes=('x0305', 'x1386', 'x1418'),
              smiles='*CCC(=O)N1CCN(Cc2sccc2C#N)CC1')

def rectifier_test():
    # name: [before, after]
    chemdex = {'phenylnaphthalene': ('c1ccc2ccccc2c1(c3ccccc3)', 'c1ccc(-c2cccc3ccccc23)cc1'),
               'benzo-azetine': ('C12CCCCC1CC2', 'C1CCC2CCCC2C1'),
               'conjoined': ('C1C2CCC2C1', 'C1CCCCC1'),
               'allene': ('C=C=C', 'C=CC'),
               'benzo-cyclopronane': ('C12CCCCC1C2', 'C1CCC2CCCC2C1'),
               'norbornane': ('C1CC2CCC1C2', 'C1CC2CCC1C2')
               }

    for name in chemdex:
        before, after = chemdex[name]
        mol = Chem.MolFromSmiles(before)
        AllChem.EmbedMolecule(mol)
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(recto.mol)
        assert gotten == after, f'{name} failed {gotten} (expected {after})'

def mixed_ring():
    mol = Chem.MolFromSmiles('c1cccc2c1CCCC2')
    mol = Chem.MolFromSmiles('C1CCCc2c1cccc2')

def cyclopentine():
    # aromatic cyclopent-ine -> cyclopentadiene
    mol = Chem.MolFromSmiles('[nH]1cccc1')
    AllChem.EmbedMolecule(mol)
    mod = Chem.RWMol(mol)
    mod.GetAtomWithIdx(0).SetAtomicNum(6)
    mol = mod.GetMol()
    recto = Rectifier(mol, atoms_in_bridge_cutoff=3).fix()
    recto.mol #C1=CCC=C1

def detriangulate_tets():
    ## test triangle
    mol = Chem.MolFromSmiles('C1CC1')
    AllChem.EmbedMolecule(mol)
    partA = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, [0, 1], addDummies=False), asMols=True)[0]
    partB = Chem.FragmentOnBonds(mol, [2], addDummies=False)

def desquarify_test():
    mol = Chem.MolFromSmiles('C1CCC1')
    AllChem.EmbedMolecule(mol)
    mol.SetProp('_Name', 'cyclobutane')
    partA = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, [0, 2], addDummies=False), asMols=True)[0]
    partA.SetProp('_Name', 'small')
    partB = Chem.FragmentOnBonds(mol, [1, 3], addDummies=False)
    partB.SetProp('_Name', 'big')


class Combine(unittest.TestCase):
    def setUp(self):
        mpro_folder = self.get_mpro_path()
        self.template = os.path.join(mpro_folder, 'template.pdb')


    def rings1(self):
        toluene = Chem.MolFromMolFile('toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        rototoluene = Chem.MolFromMolFile('rototoluene.mol')
        rototoluene.SetProp('_Name', 'rototoluene')
        v = Victor.combine(hits=[toluene, rototoluene],
                           pdb_filename=self.template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')

    def rings2(self):
        toluene = Chem.MolFromMolFile('toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneF = Chem.MolFromMolFile('transtoluene.mol')
        transtolueneF.SetProp('_Name', 'transtoluene-fuse')
        v = Victor.combine(hits=[toluene, transtolueneF],
                           pdb_filename=self.template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')

    def rings3(self):
        toluene = Chem.MolFromMolFile('toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneS = Chem.MolFromMolFile('transtoluene2.mol')
        transtolueneS.SetProp('_Name', 'transtoluene-spiro')
        # cmd.rotate('z', -90, 'rototoluene', camera=0)
        v = Victor.combine(hits=[toluene, transtolueneS],
                           pdb_filename=self.template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')


if __name__ == '__main__':
    #victor_test()
    rectifier_test()
