########################################################################################################################
__doc__ = \
    """
This is contains the class _MonsterBase to be inherited by _MonsterCommunal, then Monster
    """

########################################################################################################################

import logging
from typing import List, Optional, Dict, Sequence, Tuple, Union
from rdkit import Chem
from rdkit.Chem import rdFMCS


class _MonsterBase:
    """
    _MonsterBase -> _MonsterTracker -> _MonsterCommunal
    """

    journal = logging.getLogger('Fragmenstein')

    # overriding these seems insane.
    dummy_symbol = '*'
    dummy = Chem.MolFromSmiles(dummy_symbol)  #: The virtual atom where the targets attaches

    # settings...
    atoms_in_bridge_cutoff = 2
    # atoms_in_bridge_cutoff is how many bridge atoms can be deleted?
    # (0 = preserves norbornane, 1 = preserves adamantane)
    throw_on_discard = False
    strict_matching_mode = dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
                                 bondCompare=rdFMCS.BondCompare.CompareOrder,
                                 ringMatchesRingOnly=True,
                                 ringCompare=rdFMCS.RingCompare.StrictRingFusion,
                                 matchChiralTag=True)
    matching_modes = [
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareAny,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=False),  # this shape based matching is too permissive,
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=False),
        dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=False),
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareAny,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=True),
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=True),
        dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=True)]

    # ------------------------------------------------------------------------------------------------------------------

    def __init__(self,
                 hits: List[Chem.Mol],
                 average_position: bool=False,
                 joining_cutoff: float =5,
                 random_seed: Optional[int] = None):
        """
        Initialisation starts Monster, but it does not do any mergers or placements.
        This is changed in revision 0.6 (previously `mol` was specified for the latter)

        :param hits:
        :param average_position:
        :param joining_cutoff: joining cutoff used in "full" mode
        :param random_seed: A random seed for rdkit embedding calculations during placement

        """
        # ==== hits ===========================================
        # fix_hits: assert Chem.Mol, fix name if needed and store positions (see ``store_positions``)
        self.hits = self.fix_hits(hits)  # list of hits
        # ==== other ==========================================
        #self._debug_draw has been taken over by ``modifications`` and ``journal``
        self.average_position = average_position
        # ==== To do be filled ================================
        # -------- placement ----------------------------------
        self.initial_mol = None  # to be filled by place. The starting molecule (Chem.Mol).
        # Manually assignmnt of self.initial_mol is futile
        self.attachment = None  # place only.
        self.custom_map: Dict[str, Dict[int, int]] = {}
        # -------- common ------------------------------------
        # # ivars of type List[str]
        self.unmatched = []  # rejected hit names List[str]
        # self.matched is dynamic.  # accepted hits names List[str]
        # # ivars of type Chem.Mol or List[Chem.Mol] or Dict[Chem.Mol]
        self.modifications = {}
        self.positioned_mol = None  # final molecule
        self.joining_cutoff = joining_cutoff  # over-ridden
        self.random_seed = random_seed
        self.mol_options = []  # equally valid alternatives to self.positioned_mol
        self._collapsed_ring_offset = 0  # variable to keep track of how much to offset in ring collapse.
        # formerly:
        # self.scaffold = None  # template which may have wrong elements in place, or
        # self.chimera = None  # merger of hits but with atoms made to match the to-be-aligned mol

    def fix_hits(self, hits: List[Chem.Mol]) -> List[Chem.Mol]:
        """
        Adds the ``_Name`` Prop if needed
        asserts everything is a Chem.Mol
        calls ``store_positions``
        :param hits:
        :return:
        """
        dejavu = set()
        for hi, hit in enumerate(hits):
            if isinstance(hit, str):
                self.journal.warning(f'Hit {hi} is a string ({hit}).' +
                                     'This route is not the intended way. Trying to read it.')
                if '.mol' in hit or '.mdf' in hit:
                    hits[hi] = Chem.MolFromMolFile(hit)
                elif '.pdb' in hit:
                    hits[hi] = Chem.MolFromPDBFile(hit)
                else:
                    raise ValueError(f'Hit {hit} is not a Mol file.')
            elif isinstance(hit, Chem.Mol):
                pass
            else:
                raise ValueError(f'Hit has to be a Chem.Mol! not {type(hit)}')
            # fallback naming.
            if not hit.HasProp('_Name') or hit.GetProp('_Name').strip() == '':
                hit.SetProp('_Name', f'hit{hi}')
            elif hit.GetProp('_Name') in dejavu:
                hit.SetProp('_Name', hit.GetProp('_Name') + f'_{hi}')
            dejavu.add(hit.GetProp('_Name'))

            # ====== IMPORTANT ==========
            self.store_positions(hit)
        return hits

    def store_positions(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Saves positional data as _x, _y, _z and majorly ``_ori_i``, the original index.
        The latter gets used by ``_get_new_index``.

        :param mol:
        :return:
        """
        conf = mol.GetConformer()
        name = mol.GetProp('_Name')
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom.SetIntProp('_ori_i', i)
            atom.SetProp('_ori_name', name)
            atom.SetDoubleProp('_x', pos.x)
            atom.SetDoubleProp('_y', pos.y)
            atom.SetDoubleProp('_z', pos.z)
        return mol

    @property
    def matched(self) -> List[str]:
        """
        This is the counter to unmatched.
        It's dynamic as you never know...

        :return:
        """
        return [h.GetProp('_Name') for h in self.hits if
                h.GetProp('_Name') not in self.unmatched]

    def get_hit_by_name(self, name: str) -> Chem.Mol:
        """
        Given a name of a hit (as defined in ``_Name`` property), return the hit.
        Do note `fix_hits` will have been called, so the name may be assigned.
        :param name:
        :return:
        """
        for hit in self.hits:
            if hit.GetProp('_Name') == name:
                return hit
        else:
            raise ValueError(f"No hit with name {name}")
        

    def fix_custom_map(self,
                          custom_map: Dict[str, Union[Sequence[Tuple[int, int]], Dict[int, int]]]) \
            -> Dict[str, Dict[int, int]]:
        """
        This is duplicated in SpecialCompareAtoms,
        but will be deprecated in favour of this one.

        Make sure its Dict[str, Dict[int, int]]

        There is a bit of confusion about the custom map.
        Converts the custom map from dict of lists of 2-element tuples to dict of dicts.
        """
        if custom_map is None:
            # in Monster only not SpecalCompareAtoms!
            return {h.GetProp('_Name'): {} for h in self.hits}
        assert isinstance(custom_map, dict), 'User defined map has to be mol name to Dict[int, int]'
        for name, hit_map in custom_map.items():
            custom_map[name] = dict(hit_map)
        return custom_map

