########################################################################################################################

__doc__ = \
    """
Unmerge mapper (not inherited) but inherits `positional_mapping.GPM`
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__citation__ = ""
from ..version import __version__

########################################################################################################################


from typing import Callable, List, Dict, Tuple, Optional, Any, Sequence
from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np
import json
from .positional_mapping import GPM
from ..error import DistanceError
from collections import deque
import logging

log = logging.getLogger(__name__)


class Unmerge(GPM):
    """
    This class tries to solve the mapping problem by try all possible mappings of the target to the ligand.
    It is one of three in Monster (full merge, partial merge, unmerge.

    It is great with fragments that do not connect, but is bad when a hit has a typo.

    * the positions must overlap if any atom is mapped in two maps
    * no bond can be over 3 A

    The chosen map ``combined_map`` is a dict that goes from ``followup`` mol to ``combined`` mol which
    is the hits in a single molecule.

    Note that some molecules are discarded entirely.

    """
    max_strikes = 3  #: number of discrepancies tollerated.
    rotational_approach = True
    pick = -1  # override to pick not the lowest energy match.
    distance_cutoff = 3  #: how distance is too distant in Å

    def __init__(self,
                 followup: Chem.Mol,
                 mols: List[Chem.Mol],
                 maps: Dict[str, List[Dict[int, int]]],
                 no_discard: bool = False):
        """
        At the minute maps is a dict of hit name to a list of possible maps, wherein each map is a dict of
        atom index in the followup molecule to atom index in the target molecule,
        not the reverse as it will be soon... (see ``mcs_mapping``)

        :param followup: the molecule to place
        :type followup: Chem.Mol
        :param mols: 3D molecules
        :type mols: List[Chem.Mol]
        :param maps: can be generated outseide of Monster by ``.make_maps``
        :type maps: Dict[str, List[Dict[int, int]]]
        :param no_discard: do not allow any to be discarded
        """
        # ---- inputs ------------
        self.followup = followup
        self.no_discard = no_discard
        if self.no_discard:
            self.max_strikes = 100
        # only store non empty maps
        self.mols: List[Chem.Mol] = [mol for mol in mols if len(maps[mol.GetProp('_Name')]) != 0]
        self.maps: Dict[str, List[Dict[int, int]]] = maps
        d: int = len(mols) - len(self.mols)
        if self.no_discard and d > 0:
            raise DistanceError(hits=mols)
        # ---- to be filled ------------
        # see `.store`
        # list of indices in the followup that triggered a strike
        # debug only
        self.poisonous_indices = []
        accounted_for = set()
        self.c_map_options: List[Dict[int, int]] = []
        self.c_options: List[Chem.Mol] = []
        self.c_disregarded_options: List[List[Chem.Mol]] = []
        self.combined: Chem.Mol = Chem.Mol()
        self.combined_alternatives: List[Chem.Mol] = []
        self.combined_map: Dict[int, int] = {}
        self.disregarded: List[Chem.Mol] = [mol for mol in mols if len(maps[mol.GetProp('_Name')]) == 0]
        self.combined_bonded: Chem.Mol = Chem.Mol()
        self.combined_bonded_alternatives: List[Chem.Mol] = []
        self.combined_map_alternatives: List[Dict[int, int]] = []
        self.calculate(accounted_for)

    def calculate(self, accounted_for: set):
        """perform the calculations"""

        #  ---- sorters  --------------------
        goodness_sorter = self.goodness_sorter_factory(3)
        accounted_sorter = self.template_sorter_factory(accounted_for)
        # ---- rotate ----------------------------
        if self.rotational_approach:
            # this is the default basically
            others = deque(self.mols)
            for s in range(len(self.mols)):
                others.rotate(1)
                self.unmerge_inner(Chem.Mol(), {}, list(others), [])
        else:  # pre sort
            others = sorted(self.mols, key=accounted_sorter, reverse=True)
            self.unmerge_inner(Chem.Mol(), {}, list(others), [])
            i = sorted(range(len(self.c_options)),
                       key=goodness_sorter,
                       reverse=True)[0]
            # key changed from 0 to i subsequently: pretty sure that was a typo (even if this is an ancient route)
            for alt in self.c_disregarded_options[i]:
                aname = alt.GetProp('_Name')
                not_alt = set([o for o in others if o.GetProp('_Name') != aname])
                self.unmerge_inner(Chem.Mol(), {}, [alt] + list(not_alt), [])
        # ---- find best ------------------------------------
        if self.no_discard:
            valids = [i for i, v in enumerate(self.c_disregarded_options) if len(v) == 0]
            if len(valids) == 0:
                raise DistanceError(message='No valid mappings that do not disregard compounds.')
        else:
            valids = list(range(len(self.c_options)))
        indices = sorted(valids,
                         key=goodness_sorter,
                         reverse=True)
        i = indices[0]  # pick first for now
        if self.pick == -1:
            pass
        elif self.pick is None:  # bad option, correct it.
            self.pick = -1
        elif self.pick >= len(indices):  # override `pick` as not available
            self.pick = -1
        else:  # override applicable!
            i = indices[self.pick]
        ref = goodness_sorter(i)
        equals = [j for j in indices if goodness_sorter(j) == ref]
        if len(equals) > 1:
            log.info(f'Unmerge: There are {len(equals)} equally good mappings (this slows things down).')
        # if self._debug_draw:
        #     print(f'## Option #{i}  for combinations:')
        #     for j in range(len(self.c_options)):
        #         mol = self.c_options[j]
        #         m = self.c_map_options[j]
        #         d = self.c_disregarded_options[j]
        #         dv = self.measure_map(mol, m)
        #         print(j, [dd.GetProp('_Name') for dd in d], len(m), np.mean(dv), np.max(dv), self.offness(mol, m))
        # ----- fill ----------------------------------------------------------------
        self.combined = self.c_options[i]
        self.combined_map = self.c_map_options[i]
        self.disregarded = self.c_disregarded_options[i]
        self.combined_bonded: Chem.Mol = self.bond()
        alternative_indices = [j for j in equals if j != i]
        self.combined_alternatives = [self.c_options[j] for j in alternative_indices]
        self.combined_map_alternatives = [self.c_map_options[j] for j in alternative_indices]
        self.combined_bonded_alternatives = [self.bond(n) for n in range(len(self.combined_alternatives))]
        # ----- return ----------------------------------------------------------------
        return self

    def get_key(self, d: dict, v: Any):
        """
        Given a value and a dict and a value get the key.
        :param d:
        :param v:
        :return:
        """
        return list(d.keys())[list(d.values()).index(v)]

    @classmethod
    def make_maps(cls, target: Chem.Mol, mols: List[Chem.Mol], mode: Optional[Dict[str, Any]] = None) \
            -> Dict[str, List[Dict[int, int]]]:
        """
        This is basically if someone is using this class outside of Monster

        Returns a dictionary of key mol name and
        value a list of possible dictionary with index of an atom in target mol to the index given mol.
        Note that a bunch of mapping modes can be found in Monster init mixin class.

        :param target: the molecule to be mapped
        :param mols: the list of molecules with positional data to be mapped to
        :param mode: dict of setting for MCS step
        :return:
        """

        def get_atom_maps(molA, molB, **mode) -> List[Sequence[Tuple[int, int]]]:
            mcs = rdFMCS.FindMCS([molA, molB], **mode)
            common = Chem.MolFromSmarts(mcs.smartsString)
            matches = []
            # prevent a dummy to match a non-dummy, which can happen when the mode is super lax.
            is_dummy = lambda mol, at: mol.GetAtomWithIdx(at).GetSymbol() == '*'
            all_bar_dummy = lambda Aat, Bat: (is_dummy(molA, Aat) and is_dummy(molB, Bat)) or not (
                    is_dummy(molA, Aat) or is_dummy(molB, Bat))
            for molA_match in molA.GetSubstructMatches(common, uniquify=False):
                for molB_match in molB.GetSubstructMatches(common, uniquify=False):
                    matches.append([(molA_at, molB_at) for molA_at, molB_at in zip(molA_match, molB_match) if
                                    all_bar_dummy(molA_at, molB_at)])
            # you can map two toluenes 4 ways, but two are repeats.
            return list(set([tuple(sorted(m, key=lambda i: i[0])) for m in matches]))

        if mode is None:
            mode = dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
                        bondCompare=rdFMCS.BondCompare.CompareOrder,
                        ringMatchesRingOnly=True,
                        ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
                        matchChiralTag=True)
        maps: Dict[str, List[Dict[int, int]]] = {}
        for template in mols:
            pair_atom_maps = get_atom_maps(target, template, **mode)
            maps[template.GetProp('_Name')] = [dict(p) for p in pair_atom_maps]
        return maps

    def template_sorter_factory(self, accounted_for) -> Callable:
        """ returns the number of atoms that have not already been accounted for."""

        # key for sorting. requires outer scope ``maps`` and ``accounted_for``.
        def template_sorter(t: Chem.Mol) -> int:
            n_atoms = max([len([k for k in m if k not in accounted_for]) for m in self.maps[t.GetProp('_Name')]])
            return n_atoms

        return template_sorter

    def goodness_sorter_factory(self, offness_weight:int=3) -> Callable:
        """
        This is a factory for symmetry with template sorter... there is zero other reason for it to be so.
        """

        def goodness_sorter(i: int) -> int:
            # offness: How many bonds are too long?
            n_off_atoms: int = self.offness(self.c_options[i], self.c_map_options[i])
            return len(self.c_map_options[i]) - offness_weight * n_off_atoms

        return goodness_sorter

    def store(self, combined: Chem.Mol, combined_map: Dict[int, int], disregarded: List[Chem.Mol]):
        """
        Stores combined molecule and its map into the instance.
        """
        combined.SetProp('parts', json.dumps([m.GetProp('_Name') for m in disregarded]))
        self.c_map_options.append(combined_map)
        self.c_options.append(combined)
        self.c_disregarded_options.append(disregarded)
        return None

    def unmerge_inner(self,
                      combined: Chem.Mol,
                      combined_map: Dict[int, int],
                      others: List[Chem.Mol],
                      disregarded: List[Chem.Mol]) -> None:
        """
        Assesses a combination of maps
        rejections: unmapped (nothing maps) / unnovel (adds nothing)

        Do note that this method uses a lot of instance attributes.
        ``self.maps`` has the mapping data.
        This method combines to make ``self.combined_map``.

        :param combined:
        :param combined_map: This is passed empty the first time.
        :param others: It's a list of a deque that is rotated (if rotational_approach is True) of self.mols
        :param disregarded:
        :return:
        """
        # stop
        if len(others) == 0:
            self.store(combined=combined, combined_map=combined_map, disregarded=disregarded)
            return None
        # prevent issues.
        combined_map = dict(combined_map)
        others = list(others)
        disregarded = list(disregarded)
        # sort
        accounted_for = set(combined_map.keys())
        # parse
        other = others[0]
        oname = other.GetProp('_Name')
        ot = len(self.maps[oname])
        for oi, o_pair in enumerate(self.maps[oname]):
            o_map = dict(o_pair)
            o_present = set(o_map.keys())
            label = f'{oname} ({oi + 1}/{ot})'
            if len(o_map) == 0:
                possible_map = {}
            elif len(o_present - accounted_for) == 0:
                possible_map = {}
            elif combined.GetNumAtoms() == 0:
                possible_map = o_map
            else:
                possible_map = self.get_possible_map(other=other,
                                                     label=label,
                                                     o_map=o_map,
                                                     inter_map=self.get_positional_mapping(other, combined),
                                                     combined=combined,
                                                     combined_map=combined_map)
            # verdict
            self.judge_n_move_on(combined, combined_map, other, possible_map, others, disregarded)

    def judge_n_move_on(self, combined, combined_map, other, possible_map, others, disregarded):
        """
        The mutables need to be within their own scope

        :param combined:
        :param combined_map:
        :param other:
        :param possible_map:
        :param others:
        :param disregarded:
        :return:
        """
        if len(possible_map) == 0:
            # reject
            combined = Chem.Mol(combined)
            disregarded = [*disregarded, other]  # new obj
        else:
            # accept
            combined_map = {**combined_map, **possible_map}  # new obj
            combined = Chem.CombineMols(combined, other)  # new obj
            name = '-'.join([m.GetProp('_Name') for m in (combined, other) if m.HasProp('_Name')])
            combined.SetProp('_Name', name)
            disregarded = disregarded.copy()  # new obj
        # do inners
        accounted_for = set(combined_map.keys())
        template_sorter = self.template_sorter_factory(accounted_for)
        sorted_others = sorted(others[1:], key=template_sorter)
        self.unmerge_inner(combined, combined_map, sorted_others, disregarded)

    def get_possible_map(self,
                         other: Chem.Mol,
                         label: str,
                         o_map: Dict[int, int],  # followup -> other
                         inter_map: Dict[int, int],  # other -> combined
                         combined: Chem.Mol,
                         combined_map: Dict[int, int]) -> Dict[int, int]:
        """
        This analyses a single map (o_map) and returns a possible map

        :param other:
        :param label:
        :param o_map: followup -> other
        :param inter_map:
        :param combined:
        :param combined_map: followup -> combined
        :return: followup -> other
        """
        possible_map = {}
        strikes = 0  # x strikes is discarded
        accounted_for = set(combined_map.keys())
        for i, o in o_map.items():  # check each atom is okay
            # i = followup index
            # o = other index
            if i in accounted_for:  # this atom is accounted for. Check it is fine.
                if o in inter_map:  # this position overlaps
                    c = inter_map[o]  # equivalent index of combined
                    if c not in combined_map.values():
                        # the other atom does not contribute
                        strikes += 1
                        self.poisonous_indices.append(i)
                    elif self.get_key(combined_map, c) == i:
                        pass  # that is fine.
                    else:  # no it's a different atom
                        strikes += 1
                        self.poisonous_indices.append(i)
                else:  # this position does not overlaps. Yet atom is accounted for.
                    strikes += 1
                    self.poisonous_indices.append(i)
            elif o not in inter_map:
                # new atom that does not overlap
                possible_map[i] = combined.GetNumAtoms() + o
            elif inter_map[o] not in combined_map.values():
                # overlaps but the overlap was not counted
                possible_map[i] = combined.GetNumAtoms() + o
            else:  # mismatch!
                log.debug(f'{label} - {i} mismatch')
                strikes += 1
                self.poisonous_indices.append(i)
        if strikes >= self.max_strikes:
            return {}
        elif not self.check_possible_distances(other, possible_map, combined, combined_map,
                                               cutoff=self.distance_cutoff):
            return {}
        else:
            return possible_map

    def check_possible_distances(self, other, possible_map, combined, combined_map, cutoff=2.5):
        for i, offset_o in possible_map.items():
            unoffset_o = offset_o - combined.GetNumAtoms()
            atom = self.followup.GetAtomWithIdx(i)
            for neigh in atom.GetNeighbors():
                ni = neigh.GetIdx()
                if ni in possible_map:
                    pass  # assuming the inspiration compound was not janky
                elif ni in combined_map:
                    separation = self.get_inter_distance(other, combined, unoffset_o, combined_map[ni])
                    # removing the distance cutoff is a bad idea as the next atoms along with be stretched
                    # but on others it is fine.
                    if separation > cutoff * 2:  # very bad.
                        return False
                    elif separation > cutoff:  # moderately bad.
                        self.poisonous_indices.append(i)
                        self.poisonous_indices.append(ni)
                        return False
                else:
                    pass  # unmapped neighbor
        return True

    def bond(self, idx: Optional[int] = None) -> Chem.Mol:
        """
        Add bonds. As in the verb 'to bond' ...
        """
        if idx is None:
            # calculating  self.combined_bonded
            mol = self.combined
            mapping = self.combined_map
        else:
            # calculating one of self.combined_bonded_alternatives
            mol = self.combined_alternatives[idx]
            mapping = self.combined_map_alternatives[idx]
        putty = Chem.RWMol(mol)
        for fi, ci in mapping.items():
            fatom = self.followup.GetAtomWithIdx(fi)
            for neigh in fatom.GetNeighbors():
                ni = neigh.GetIdx()
                if ni not in mapping:
                    continue
                nci = mapping[ni]
                bond_type = self.followup.GetBondBetweenAtoms(fi, ni).GetBondType()
                if not putty.GetBondBetweenAtoms(ci, nci):
                    putty.AddBond(ci, nci, bond_type)
                else:
                    putty.GetBondBetweenAtoms(ci, nci).SetBondType(bond_type)
        return putty.GetMol()

    def get_inter_distance(self, molA: Chem.Mol, molB: Chem.Mol, idxA: int, idxB: int) -> float:
        def get_pos(mol, idx):
            conf = mol.GetConformer()
            return np.array(conf.GetAtomPosition(idx))

        return np.linalg.norm(get_pos(molA, idxA) - get_pos(molB, idxB))

    def measure_map(self, mol: Chem.Mol, mapping: Dict[int, int]) -> np.array:
        """
        Returns a vector with the distances but not of length len(mapping)
        This used by offness to score how bad the mapping is

        :param mol:
        :param mapping: followup to comined
        :return:
        """
        conf = mol.GetConformer()
        atomic_distances = np.array([])
        for followup_atom_idx, template_atom_idx in mapping.items():
            # followup mol is `victor.monster.initial_mol`
            followup_atom = self.followup.GetAtomWithIdx(followup_atom_idx)
            for followup_neigh_atom in followup_atom.GetNeighbors():
                followup_neigh_idx = followup_neigh_atom.GetIdx()
                if followup_neigh_idx not in mapping:
                    continue
                template_neigh_idx = mapping[followup_neigh_idx]
                # calculate euclidean distance between the two atoms in the followup mol
                a = np.array(conf.GetAtomPosition(template_atom_idx))
                b = np.array(conf.GetAtomPosition(template_neigh_idx))
                atomic_distances = np.append(atomic_distances, np.linalg.norm(a - b))
        return atomic_distances

    def offness(self, mol: Chem.Mol,
                mapping: Dict[int, int],
                cutoff_distance: float = 2.5) -> int:
        """
        How many bonds are too long?

        :param mol:
        :param mapping:
        :return:
        """
        d = self.measure_map(mol, mapping)
        # return np.linalg.norm(d - 1.5)/(d.size*0.5) # 1.5 ang
        return sum(d > cutoff_distance)
