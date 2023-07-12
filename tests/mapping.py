import unittest
from rdkit import Chem
from fragmenstein import Monster
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from fragmenstein.monster.mcs_mapping import SpecialCompareAtoms
# ======================================================================================================================

class Mappings(unittest.TestCase):
    """
    These test monster.mcs_mapping
    """
    toluene = Chem.MolFromSmiles('Cc1ccccc1')
    toluene.SetProp('_Name', 'toluene')

    benzyl = Chem.MolFromSmiles('*c1ccccc1')
    benzyl.SetProp('_Name', 'benzyl')

    methylpyrylium = Chem.MolFromSmiles('Cc1c[o+]ccc1')
    methylpyrylium.SetProp('_Name', 'methylpyrylium')
    # Monster.draw_nicely(None, methylpyrylium)

    methylpyridine = Chem.MolFromSmiles('Cc1ncccc1')
    methylpyridine.SetProp('_Name', 'methylpyridine')
    # Monster.draw_nicely(None, methylpyridine)


    def test_not_to_dummy(self):
        from rdkit.Chem import rdFMCS
        params = rdFMCS.MCSParameters()  # noqa
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms()
        compare = [self.benzyl, self.toluene]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        print(res.smartsString)
        self.assertEqual(res.numAtoms, 6)  # there are 7 atoms, but only 6 are mapped as the dummy is excluded

    def test_user_map(self):
        """
        Map the azo with the oxo in pyridine and pyrylium.
        Using the methyl versions, the methyl does not map.
        """
        from rdkit.Chem import rdFMCS
        custom_map = {'methylpyridine': {2: 3}}
        params = rdFMCS.MCSParameters()
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms(custom_map=custom_map)
        compare = [self.methylpyridine, self.methylpyrylium]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        self.assertEqual(res.numAtoms, 6, 'Mapping  the azo with the oxo means no methyl res:' +
                                         f' 6 matched, not {res.numAtoms}')
        full_maps = params.AtomTyper.get_valid_matches(params.AtomCompareParameters,
                                                      Chem.MolFromSmarts(res.smartsString),
                                                      self.methylpyridine, self.methylpyrylium
                                                      )
        # this is wrong: you can map the methyls to either side (ortho) of the heteroatom
        # self.assertEqual(len(full_maps), 1, f'There is only one way to map it not {len(full_maps)} ({full_maps})')
        for full_map in full_maps:
            for h, f in custom_map['methylpyridine'].items():
                self.assertEqual(dict(full_map)[h], f, f'hit atom {h} => {f} not {dict(full_map)[h]}')

    def test_user_negmap(self):
        """
        Ban index 1 on the azo=oxo mapping
        """
        from rdkit.Chem import rdFMCS
        custom_map = {'methylpyridine': {2: 3, 1: -2}}
        params = rdFMCS.MCSParameters()
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms(custom_map=custom_map)
        compare = [self.methylpyridine, self.methylpyrylium]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        self.assertEqual(res.numAtoms, 5, 'ought to be azo=oxo and no index 1')
        full_maps = params.AtomTyper.get_valid_matches(params.AtomCompareParameters,
                                                       Chem.MolFromSmarts(res.smartsString),
                                                       self.methylpyridine, self.methylpyrylium
                                                       )
        # this is wrong: you can map the methyls to either side (ortho) of the heteroatom
        # self.assertEqual(len(full_maps), 1, f'There is only one way to map it not {len(full_maps)} ({full_maps})')
        for full_map in full_maps:
            self.assertNotIn(1, dict(full_map).keys())

    def test_flipper(self):
        from fragmenstein.monster.mcs_mapping import flip_mapping
        # Sequence[Tuple[int, int]]
        flipped = flip_mapping([(1, 2)])
        self.assertIsInstance(flipped, (list, tuple))


    def test_user_negmap2(self):
        """
        Assign followup index 1 to another molecule ('nullium' as its nothing) in the azo=oxo mapping
        thus preventing its mapping to methylpyridine
        """
        from rdkit.Chem import rdFMCS
        custom_map = {'methylpyridine': {2: 3}, 'nullium': {1: 1}}
        params = rdFMCS.MCSParameters()
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms(custom_map=custom_map)
        compare = [self.methylpyridine, self.methylpyrylium]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        self.assertEqual(res.numAtoms, 5, 'ought to be azo=oxo and no index 1')
        full_maps = params.AtomTyper.get_valid_matches(params.AtomCompareParameters,
                                                       Chem.MolFromSmarts(res.smartsString),
                                                       self.methylpyridine, self.methylpyrylium
                                                       )
        # this is wrong: you can map the methyls to either side (ortho) of the heteroatom
        # self.assertEqual(len(full_maps), 1, f'There is only one way to map it not {len(full_maps)} ({full_maps})')
        for full_map in full_maps:
            self.assertNotIn(1, dict(full_map).values())

    def test_issue42(self):
        RDLogger.DisableLog('rdApp.warning')  # shut up about the hydrogens already
        # hydrogenated indole w/ a methyl on the r6
        hit = Chem.MolFromSmiles('N1CCC2C1CC(C)CC2')
        AllChem.EmbedMolecule(hit)
        hit.SetProp('_Name', 'foo')
        monstah = Monster([hit])
        # hydrogenated benzothiazole w/ no methyl
        monstah.place_smiles('N1CSC2C1CCCC2')
        self.assertEqual('foo.0+foo.1+foo.2+foo.3+foo.4+foo.5+foo.6+foo.8+foo.9',
                         '+'.join([o[0] if o else 'NA' for o in monstah.origin_from_mol()]))



if __name__ == '__main__':
    unittest.main()
