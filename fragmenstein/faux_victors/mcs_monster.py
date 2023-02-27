import json, logging
from rdkit import Chem, Geometry
from rdkit.Chem import rdFMCS, AllChem  # noqa
from typing import Dict
from molecular_rectifier import Rectifier
from ..monster import Monster


class MCSMerger(Monster):
    """
    see ``combine()`` for ops.
    """
    journal = logging.getLogger()

    extreme_embed_args = [('maxAttempts', 0, 5),
                          ('numZeroFail', 1, 3),
                          ('forceTol', 0.001, 1),
                          ('ignoreSmoothingFailures', False, True),  # max([True, False]) == True
                          ]

    def __init__(self, hits: Chem.Mol, mcs_parameters=rdFMCS.MCSParameters()):
        self.hits = hits
        assert len(hits) == 2, 'Only two hits requirement, not more, not less'
        self.fore = hits[0]
        self.aft = hits[1]
        self.mcs_parameters = mcs_parameters
        self.common = self.get_common()
        self.aft2fore_mapping = self.get_aft2fore_mapping()
        # added in combine
        self.mod = Chem.RWMol
        self.flat_combo = Chem.Mol
        self.rectified = Chem.Mol
        self.positioned_mol = Chem.Mol
        # these are to match Monster
        self.minimized_mol = Chem.Mol
        self.mol_options = []
        self.unmatched = []

    def get_common(self) -> Chem.Mol:
        res = rdFMCS.FindMCS([self.fore, self.aft], self.mcs_parameters)
        return Chem.MolFromSmarts(res.smartsString)

    def get_aft2fore_mapping(self) -> Dict[int, int]:
        return dict(zip(self.aft.GetSubstructMatch(self.common), self.fore.GetSubstructMatch(self.common)))

    def combine(self, embed_mode: int = 0):
        """
        Combine fore and aft into a single molecule (see below).
        Will add property 'provenance' to atoms, with values 'fore' or 'aft' or 'common'.
        Will add property 'provenance' to bonds, with values 'fore', 'aft', 'common' or 'semicommon' (different bond order).
        This property is different from `_Origin` / `_Origins` and `_BondProvenance`...

        The argument ``embed_mode`` (int for now) controns if to use
        0 = ``AllChem.EmbedMolecule`` or
        1 = ``AllChem.ConstrainedEmbed``

        * ``.mod`` is the combined molecule (RWMol).
        * ``.flat_combo`` is ``.mod`` without Hs and coordinates.
        * ``.rectified`` is ``.flat_combo`` after rectification.
        * ``.positioned_mol`` is ``.rectified`` with Hs and embedded coordinates.
        """
        self.mark(self.fore, 'fore_')
        self.mod = Chem.RWMol(self.fore)
        # copy atoms and bonds from aft to fore-copy
        self._add_atoms()
        # copy bonds from aft to fore-copy
        self._add_bonds()
        # now mod > flat_combo > rectified > positioned_mol
        # clean up
        self.positioned_mol = self._clean_up()
        # embed
        if embed_mode == 0:
            self._embed()
        else:
            self._con_embed()
        # self._con_embed()
        self.positioned_mol.SetProp('_Name', self.merged_name)
        self.convert_to_origins(self.positioned_mol)
        return self.positioned_mol

    @property
    def merged_name(self):
        if not self.fore.HasProp('_Name') and not self.aft.HasProp('_Name'):
            return 'merger'
        if not self.fore.HasProp('_Name'):
            self.fore.SetProp('First')
        if not self.aft.HasProp('_Name'):
            self.aft.SetProp('Second')
        return self.fore.GetProp('_Name') + self.aft.GetProp('_Name')

    def mark(self, mol: Chem.Mol, prefix=''):
        """
        Adds the fore_i and aft_i properties to atoms.
        """
        for atom in mol.GetAtoms():
            atom.SetIntProp(prefix + 'i', atom.GetIdx())

    def _add_atoms(self):
        """
        Add atoms from aft to mod, a copy of fore.
        Called by combine().
        Will add property 'provenance' to atoms, with values 'fore', 'aft' or 'common'.
        If the atom is a product of rectification then it would be `rectified`.

        :return: None
        """
        atom: Chem.Atom
        # ### mark originals
        for atom in self.mod.GetAtoms():
            atom.SetProp('provenance', 'fore')
        # ### copy atoms
        self.mod.BeginBatchEdit()
        for atom in self.aft.GetAtoms():
            i = atom.GetIdx()
            if i in self.aft2fore_mapping:
                self.mod.GetAtomWithIdx(self.aft2fore_mapping[i]).SetProp('provenance', 'common')
                continue
            if atom.GetAtomicNum() == 1:
                # skip H
                # zahl of dummy is 0: these are left - wise?
                continue
            new_i = self.mod.AddAtom(Chem.Atom(atom.GetAtomicNum()))
            new_atom = self.mod.GetAtomWithIdx(new_i)
            new_atom.SetProp('provenance', 'aft')
            new_atom.SetIntProp('aft_i', i)
            self.aft2fore_mapping[i] = new_i
        self.mod.CommitBatchEdit()

    def _add_bonds(self):
        """
        Add bonds from aft to mod, a copy of fore.
        Called by combine().
        Will add property 'provenance' to bonds, with values 'fore' or 'aft' or 'common' or 'semicommon'.
        The latter for bonds that are common but with different bond order.
        If the bond is a product of rectification then it would be `rectified`.

        :return: None
        """
        bond: Chem.Bond
        for bond in self.mod.GetBonds():
            bond.SetProp('provenance', 'fore')
            bond.SetIntProp('fore_i', bond.GetIdx())
        # ### copy bonds
        for aft_bond in self.aft.GetBonds():
            aft_begin_idx = aft_bond.GetBeginAtomIdx()
            if self.aft.GetAtomWithIdx(aft_begin_idx).GetAtomicNum() == 1:
                # skip H
                continue
            aft_end_idx = aft_bond.GetEndAtomIdx()
            if self.aft.GetAtomWithIdx(aft_end_idx).GetAtomicNum() == 1:
                # skip H
                continue
            fore_begin_idx = self.aft2fore_mapping[aft_begin_idx]
            fore_end_idx = self.aft2fore_mapping[aft_end_idx]
            mod_bond = self.mod.GetBondBetweenAtoms(fore_begin_idx, fore_end_idx)  # noqa
            if mod_bond is None:
                # the bond did not exist
                total_bonds: int = self.mod.AddBond(fore_begin_idx, fore_end_idx, aft_bond.GetBondType())  # noqa
                mod_bond = self.mod.GetBondBetweenAtoms(fore_begin_idx, fore_end_idx)  # noqa
                mod_bond.SetProp('provenance', 'aft')
            elif aft_bond.GetBondTypeAsDouble() > mod_bond.GetBondTypeAsDouble():
                # the bond existed but wasn't as reduced
                mod_bond.SetBondType(aft_bond.GetBondType())
                mod_bond.SetProp('provenance', 'semicommon')
            else:
                mod_bond.SetProp('provenance', 'common')

    def _clean_up(self):
        """
        `.flat_combo`` is ``.mod`` withoout Hs and coordinates.
        `.rectified`` is ``.flat_combo`` after rectification.
        the returned will be `.positioned_mol`` ie. ``.rectified`` with Hs.
        This will gain coordinates in `.embed()`.

        :return:
        """
        combo = self.mod.GetMol()
        self.remove_chiral_tags(combo)
        self.flat_combo = AllChem.RemoveHs(combo, sanitize=False, updateExplicitCount=True)
        AllChem.Compute2DCoords(self.flat_combo)  # why is this needed? Is rectifier 3D independent?
        self.rectified = Rectifier(self.flat_combo).fix().mol
        # get only first largest fragment (in case it's split)
        largest_frag = sorted(AllChem.GetMolFrags(AllChem.RemoveHs(self.rectified), asMols=True),
                              key=lambda mol: mol.GetNumHeavyAtoms(),
                              reverse=True)[0]
        hydroed = AllChem.AddHs(largest_frag)
        AllChem.EmbedMolecule(hydroed)
        for atom in hydroed.GetAtoms():
            if not atom.HasProp('provenance'):
                atom.SetProp('provenance', 'rectified')
        for bond in hydroed.GetBonds():
            if not bond.HasProp('provenance'):
                bond.SetProp('provenance', 'rectified')
        hydroed.SetProp('_Name', self.merged_name)
        self.remove_chiral_tags(hydroed)
        return hydroed  # will be saved as self.positioned_mol

    def remove_chiral_tags(self, mol: Chem.Mol):
        for atom in mol.GetAtoms():
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

    def is_valid(self, atom: Chem.Atom, trial: int):
        """
        Is this a constrainable atom?
        """
        if atom.GetAtomicNum() == 1:
            return False
        if not atom.HasProp('fore_i'):
            return False
        if trial < 1:
            return True
        for bond in atom.GetBonds():
            if bond.GetProp('provenance') == ('semicommon'):
                return False
        return trial == 1 or atom.GetProp('provenance') == 'common'

    def _embed(self, trial=0, **embed_args) -> None:
        """
        At trial 0, all fore atoms are okay
        At trial 1, only fore atoms with bonds with original bonding are okay
        At trial 2, only common atoms with bonds with original bonding are okay
        """
        self.journal.debug(f'Embedding trial {trial}')
        combo2fore = {atom.GetIdx(): atom.GetIntProp('fore_i') for atom in self.positioned_mol.GetAtoms() if
                      self.is_valid(atom, trial)}
        fore_conf = self.fore.GetConformer()
        coord_map: Dict[int, Geometry.Point3D] = {ci: fore_conf.GetAtomPosition(fi) for ci, fi in combo2fore.items()}
        conf_i = AllChem.EmbedMolecule(self.positioned_mol,
                                       **{**dict(clearConfs=True,
                                                 coordMap=coord_map,
                                                 ignoreSmoothingFailures=True),
                                          **embed_args
                                          }
                                       )

        if conf_i != -1:
            # add selection & return
            self.positioned_mol.__sssAtoms = list(coord_map)
            for atom in self.positioned_mol.GetAtoms():
                i: int = atom.GetIdx()
                atom.SetBoolProp('constrained', i in coord_map)
            self.positioned_mol.SetIntProp('embedding_trial', trial)
            # align
            rmsd: float = AllChem.AlignMol(self.positioned_mol,
                                           self.fore,
                                           atomMap=list(combo2fore.items())
                                           )
            # print('.... ', list(combo2fore.items()), rmsd)
            return
        # repeat...
        if trial == 2:
            # getting with rdkit default
            for key, default, new in self.extreme_embed_args:
                embed_args[key] = max(embed_args.get(key, default), new)
        if trial == 3:
            raise ValueError('Embedding failed')
        self._embed(trial + 1, **embed_args)

    def _con_embed(self, trial: int = 0) -> None:
        slimmed = Chem.RWMol(self.positioned_mol)
        keepers = []
        unkeepers = []
        for atom in slimmed.GetAtoms():
            if self.is_valid(atom, trial):
                keepers.append(atom.GetIdx())
            else:
                unkeepers.append(atom.GetIdx())

        slimmed.BeginBatchEdit()
        for i in unkeepers:
            slimmed.RemoveAtom(i)
        slimmed.CommitBatchEdit()
        fore_conf = self.fore.GetConformer()
        conf = Chem.Conformer(slimmed.GetNumAtoms())
        for atom in slimmed.GetAtoms():
            point: Geometry.Point3D = fore_conf.GetAtomPosition(atom.GetIntProp('fore_i'))
            conf.SetAtomPosition(atom.GetIdx(), point)
        slimmed.AddConformer(conf)
        try:
            AllChem.ConstrainedEmbed(self.positioned_mol, slimmed, useTethers=True, clearConfs=True)
            self.positioned_mol.__sssAtoms = keepers
            # align
            rmsd: float = AllChem.AlignMol(self.positioned_mol,
                                           self.fore,
                                           atomMap=[(ci, self.positioned_mol.GetAtomWithIdx(ci).GetIntProp('fore_i'),)
                                                    for ci in keepers]
                                           )
            for atom in self.positioned_mol.GetAtoms():
                i: int = atom.GetIdx()
                atom.SetBoolProp('constrained', i in keepers)
            self.positioned_mol.SetIntProp('embedding_trial', trial)
            return
        except ValueError:
            if trial == 2:
                raise ValueError('Embedding failed')
            self._con_embed(trial + 1)

    def convert_to_origins(self, mol: Chem.Mol, both=False):
        """
        In this class the 'provenance' are marked as common, aft, fore, rectified (semicommon is bond)
        """
        fore_name = self.fore.GetProp('_Name')
        aft_name = self.fore.GetProp('_Name')
        for atom in mol.GetAtoms():
            origins = []
            if atom.HasProp('fore_i'):
                origins.append(fore_name + '.' + atom.GetProp('fore_i'))
            if both and atom.HasProp('aft_i'):
                origins.append(aft_name + '.' + atom.GetProp('aft_i'))
            atom.SetProp('_Origin', json.dumps(origins))
