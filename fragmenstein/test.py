from . import Fragmenstein
from rdkit import Chem

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
                  hitnames=('x0434','x0540'))

def nasty2_test():
    test_molecule(name='BEN-VAN-c98-4',
                  smiles='*C(=N)CN1CCN(Cc2ccc(-c3cc(CC)ncn3)c(F)c2)CC1',
                  hitnames='x0692,x0770,x0995'.split(','))