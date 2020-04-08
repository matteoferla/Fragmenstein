# Fragmenstein
Scaffold hopping between bound compounds by stitching them together like a reanimated corpse.

To be filled out!

[https://github.com/matteoferla/SARS-CoV-2_CL3_covalent_docking](my messy code for Covid19 moonshot).

    hits = [Chem.MolFromMolFile(f'../Mpro/Mpro-{i}_0/Mpro-{i}_0.mol') for i in ('x0692', 'x0305', 'x1249')]
    probe = Chem.MolFromSmiles('CCNc1nc(CCS)c(C#N)cc1CN1C(CCS)CN(C(C)=O)CC1')
    monster = Fragmenstein(probe, hits)
    #monster = Fragmenstein(probe, hits, draw=True) for verbosity in a Jupyter notebook
    monster.make_pse('test.pse')
    
    display(monster.scaffold)
    display(monster.chimera) # merger of hits but with atoms made to match the to-be-aligned mol
    display(monster.positioned_mol) # followup aligned
    
    # further alignments... badly written way of doing this.
    monster.initial_mol = new_mol
    aligned = monster.place_probe(new_mol)