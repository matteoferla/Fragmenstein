# these are not tests. But cases as a human I would not know how to solve.

def ring_mess():
    """
    This is a case of two rings that are close and result in a weird molecule due to overzealous absorption.
    """
    hits = [Mac1.get_mol('diamond-x0711_A'), Mac1.get_mol('diamond-x0689_A')]
    Victor.monster_throw_on_discard = True
    Victor.quick_reanimation = True
    victor = Victor(hits=hits,
                    pdb_block=Mac1.get_template(),
                    )
    victor.combine()
    mol = victor.monster.modifications['expanded']

    from molecular_rectifier import Rectifier

    recto = Rectifier(mol)
    recto.fix()
    # recto.modifications[3] is good as gold
    # recto.modifications[4] is a mess