from fragmenstein import Fragmenstein, Victor, Egor

from rdkit import Chem
import logging


def test_molecule(name, smiles, hitnames):
    hits = [Chem.MolFromMolFile(f'../Mpro/Mpro-{i}_0/Mpro-{i}_0.mol') for i in hitnames]
    followup = Chem.MolFromSmiles(smiles)
    r = Chem.MolFromMolFile(f'../Mpro/Mpro-{hitnames[0]}_0/Mpro-{hitnames[0]}_0_SG.mol')
    f = Fragmenstein(followup, hits, attachment=r)
    f.make_pse(f'test_{name}.pse')
    print(f.logbook)


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


def victor_test():
    mpro_folder = '/Users/matteo/Coding/rosettaOps/Mpro'
    hits = [Chem.MolFromMolFile(f'{mpro_folder}/Mpro-{i}_0/Mpro-{i}_0.mol') for i in ('x0305', 'x1386', 'x1418')]
    Victor.enable_stdout(logging.DEBUG)
    Victor(smiles='*CCC(=O)N1CCN(Cc2sccc2C#N)CC1',
           hits=hits,
           pdb_filename=f'{mpro_folder}/Mpro-x1386_0/Mpro-x1386_0.mol',
           long_name='DAV-CRI-d1e-2_ACR',
           ligand_resn='LIG',
           ligand_resi='1B',
           covalent_resn='CYS', covalent_resi='145'
           )


if __name__ == '__main__':
    victor_test()
