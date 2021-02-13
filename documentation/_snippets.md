These are snippets that need to go somewhere
# 1

To make Victor carp at any exception, do

    Victor.error_to_catch = ()
    
This is because `try`/`except` can accept multiple error classes as a tuple.

# 2

`ConnectionError` is meant for comms connection errors, but here it is raised as bad bonding.

# 3

Show the steps

    from IPython.display import display
    
    monster = Monster([fore, aft]).merge()
    for m in monster.modifications:
        display(m)
        
Show multiple RDKit mols in nglview

    import nglview as nv
    from io import StringIO

    def show_mols(*mols: Chem.Mol) -> nv.NGLWidget:
        view = nv.NGLWidget()
        for mol in mols:
            fh = StringIO(Chem.MolToPDBBlock(mol))
            view.add_component(fh, ext='pdb')
        return view

# X

One issue is that the ring mergers are not transitive.
Hit order results in different cases:

    ortho = Chem.MolFromSmiles('Cc1cc(O)ccc1')
    meta = Chem.MolFromSmiles('Cc1ccc(O)cc1')
    AllChem.EmbedMolecule(ortho)
    AllChem.EmbedMolecule(meta)
    Chem.rdMolAlign.AlignMol(prbMol=ortho,
                             refMol=meta,
                             atomMap=[(0,0), (1,1)],
                             weights=[1,1],
                            )
    show_mols(ortho, meta)
    
These two cresol rings are perpendicular. The merger produces the same compound, methylcatecol
but one has a more distorted conformation.

    Monster([ortho, meta]).merge().positioned_mol  # decent
    Monster([meta, ortho]).merge().positioned_mol  # indecent
    
What are the implications for this?