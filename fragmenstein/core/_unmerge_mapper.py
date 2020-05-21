from typing import Callable, List, Dict, Tuple
from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np

from ._positional_mapping import GPM
from collections import deque

class Unmerge(GPM):
    """
    This class tries to solve the mapping problem by try all possible mappings of the target to the ligand.
    It is one of three in Fragmenstein (full merge, partial merge (preferred), unmerge.

    It is great with fragments that do not connect, but is bad when a hit has a typo.

    * the positions must overlap if any atom is mapped in two maps
    * no bond can be over 3 A

    The chosen map ``combined_map`` is a dict that goes from ``followup`` mol to ``combined`` mol which
    is the hits in a single molecule.

    Note that some molecules are discarded entirely.

    """
    max_strikes = 3  #: number of discrepancies tollerated.
    rotational_approach = True

    def __init__(self, followup: Chem.Mol, mols: List[Chem.Mol], maps: Dict[str, List[Dict[int, int]]],
                 _debug_draw: bool = False):
        """


        :param followup: the molecule to place
        :type followup: Chem.Mol
        :param mols: 3D molecules
        :type mols: List[Chem.Mol]
        :param maps: can be generated outseide of Fragmenstein by ``.make_maps``.
        :type maps: Dict[List[Dict[int, int]]]
        :param _debug_draw:
        """
        self.followup = followup
        self.mols = mols
        self.maps = maps
        self._debug_draw = _debug_draw
        accounted_for = set()
        self.c_map_options = []
        self.c_options = []
        self.c_disregarded_options = []
        # sorters
        goodness_sorter = lambda i: len(self.c_map_options[i]) - self.offness(self.c_options[i], self.c_map_options[i])
        accounted_sorter = self.template_sorter_factory(accounted_for)
        # rotate
        if self.rotational_approach:
            others = deque(self.mols)
            for s in range(len(self.mols)):
                others.rotate(1)
                if self._debug_draw:
                    print(f"Rotated, new first : {others[0].GetProp('_Name')}")
                self.unmerge_inner(Chem.Mol(), {}, list(others), [])
        else: # pre sort
            others = sorted(self.mols, key=accounted_sorter, reverse=True)
            self.unmerge_inner(Chem.Mol(), {}, list(others), [])
            i = sorted(range(len(self.c_options)),
                             key=goodness_sorter,
                             reverse=True)[0]
            for alt in self.c_disregarded_options[0]:
                aname = alt.GetProp('_Name')
                not_alt = set([o for o in others if o.GetProp('_Name') != aname])
                self.unmerge_inner(Chem.Mol(), {}, [alt]+list(not_alt), [])
        # find best
        indices = sorted(range(len(self.c_options)),
                         key=goodness_sorter,
                         reverse=True)
        i = indices[0]
        if self._debug_draw or 1==1:
            print(f'## Option #{i}  for combinations:')
            for j in range(len(self.c_options)):
                mol = self.c_options[j]
                m = self.c_map_options[j]
                d = self.c_disregarded_options[j]
                dv = self.measure_map(mol, m)
                print(j, [dd.GetProp('_Name') for dd in d], len(m), np.mean(dv), np.max(dv), self.offness(mol, m))
        self.combined = self.c_options[i]
        self.combined_map = self.c_map_options[i]
        self.disregarded = self.c_disregarded_options[i]
        self.combined_bonded = self.bond()

    get_key = lambda self, d, v: list(d.keys())[list(d.values()).index(v)]

    @classmethod
    def make_maps(cls, target: Chem.Mol, mols: List[Chem.Mol]) -> Dict[str, List[Dict[int, int]]]:
        """
        This is basically if someone is using this class outside of Fragmenstein

        Returns a dictionary of key mol name and
        value a list of possible dictionary with idex of target to the index given mol.
        Note that a bunch of mapping modes can be found in Fragmenstein init mixin class.

        :param target: the molecule to be mapped
        :param mols: the list of molecules with positional data to be mapped to
        :return:
        """
        def get_atom_maps(molA, molB, **mode) -> List[List[Tuple[int, int]]]:
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
            return list[set([tuple(sorted(m, key=lambda i: i[0])) for m in matches])]

        maps = {}
        for template in mols:
            pair_atom_maps = get_atom_maps(target, template,
                                                 atomCompare=rdFMCS.AtomCompare.CompareElements,
                                                 bondCompare=rdFMCS.BondCompare.CompareOrder,
                                                 ringMatchesRingOnly=True,
                                                 ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
                                                 matchChiralTag=True)
            maps[template.GetProp('_Name')] = [dict(p) for p in pair_atom_maps]
        return maps

    def template_sorter_factory(self, accounted_for) -> Callable:
        """ returns the number of atoms that have not already been accounted for."""
        # key for sorting. requires outer scope ``maps`` and ``accounted_for``.
        def template_sorter(t: Chem.Mol) -> int:
            n_atoms = max([len([k for k in m if k not in accounted_for]) for m in self.maps[t.GetProp('_Name')]])
            return n_atoms

        return template_sorter

    def store(self, combined: Chem.Mol, combined_map: Dict[int, int], disregarded: List[Chem.Mol]):
        self.c_map_options.append(combined_map)
        self.c_options.append(combined)
        self.c_disregarded_options.append(disregarded)
        return None

    def unmerge_inner(self,
                      combined: Chem.Mol,
                      combined_map: Dict[int, int],
                      others: List[Chem.Mol],
                      disregarded: List[Chem.Mol]) -> None:
        # stop
        if len(others) == 0:
            if self._debug_draw:
                print('************************ (stored)')
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
                if self._debug_draw:
                    print(f'{label} unmapped')
                possible_map = {}
            elif len(o_present - accounted_for) == 0:
                if self._debug_draw:
                    print(o_present, accounted_for)
                    print(f'{label} unnovel')
                possible_map = {}
            elif combined.GetNumAtoms() == 0:
                if self._debug_draw:
                    print(f'{label} first one')
                possible_map = o_map
            else:
                if self._debug_draw:
                    print(f'{label} assessment')
                possible_map = self.get_possible_map(other=other,
                                                     label=label,
                                                     o_map=o_map,
                                                     inter_map=self.get_positional_mapping(other, combined),
                                                     combined=combined,
                                                     combined_map=combined_map)
            # verdict
            self.judge_n_move_on(possible_map)

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
            if self._debug_draw:
                print('>> reject')
            disregarded = disregarded + [other]  # new obj
        else:
            # accept
            if self._debug_draw:
                print(f'>> accept: {possible_map}')
            combined_map = {**combined_map, **possible_map} # new obj
            combined = Chem.CombineMols(combined, other) # new obj
            name = '-'.join([m.GetProp('_Name') for m in (combined, other) if m.HasProp('_Name')])
            combined.SetProp('_Name', name)
            disregarded = list(disregarded) # new obj
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
                        if self._debug_draw:
                            print(f'{label} - {i} accounted, but no contrib')
                        strikes += 1
                    elif self.get_key(combined_map, c) == i:
                        pass  # that is fine.
                    else:  # no it's a different atom
                        if self._debug_draw:
                            print(f'{label} - {i} accounted, diff atom')
                        strikes += 1
                else:  # this position does not overlaps. Yet atom is accounted for.
                    if self._debug_draw:
                        print(f'{label} - {i} accounted, no overlap')
                    strikes += 1
            elif o not in inter_map:
                # new atom that does not overlap
                possible_map[i] = combined.GetNumAtoms() + o
            elif inter_map[o] not in combined_map.values():
                # overlaps but the overlap was not counted
                possible_map[i] = combined.GetNumAtoms() + o
            else:  # mismatch!
                if self._debug_draw:
                    print(f'{label} - {i} mismatch')
                strikes += 1
        if strikes >= self.max_strikes:
            if self._debug_draw:
                print(f'{label} got {strikes} strikes')
            return {}
        elif not self.check_possible_distances(other, possible_map, combined, combined_map):
                if self._debug_draw:
                    print(f'{label} gives too long bonds')
                return {}
        else:
            return possible_map

    def check_possible_distances(self, other, possible_map, combined, combined_map, cutoff=3):
        for i, offset_o in possible_map.items():
            unoffset_o = offset_o - combined.GetNumAtoms()
            atom = self.followup.GetAtomWithIdx(i)
            for neigh in atom.GetNeighbors():
                ni = neigh.GetIdx()
                if ni in possible_map:
                    pass # assuming the inspiration compound was not janky
                elif ni in combined_map:
                    if self.get_inter_distance(other, combined, unoffset_o, combined_map[ni]) > cutoff:
                        return False
                else:
                    pass # unmapped neighbor
        return True

    def bond(self):
        putty = Chem.RWMol(self.combined)
        for fi, ci in self.combined_map.items():
            fatom = self.followup.GetAtomWithIdx(fi)
            for neigh in fatom.GetNeighbors():
                ni = neigh.GetIdx()
                if ni not in self.combined_map:
                    continue
                nci = self.combined_map[ni]
                bond_type = self.followup.GetBondBetweenAtoms(fi, ni).GetBondType()
                if not putty.GetBondBetweenAtoms(ci, nci):
                    if self._debug_draw:
                        print(fi, ni, 'bond new')
                    putty.AddBond(ci, nci, bond_type)
                else:
                    if self._debug_draw:
                        print(fi, ni, 'bond added')
                    putty.GetBondBetweenAtoms(ci, nci).SetBondType(bond_type)
        return putty.GetMol()

    def get_inter_distance(self, molA:Chem.Mol, molB: Chem.Mol, idxA: int, idxB: int) -> np.float:
        def get_pos(mol, idx):
            conf = mol.GetConformer()
            return np.array(conf.GetAtomPosition(idx))
        return np.linalg.norm(get_pos(molA, idxA) - get_pos(molB, idxB))

    def measure_map(self, mol: Chem.Mol, mapping: Dict[int, int]) -> np.array:
        """

        :param mol:
        :param mapping: followup to comined
        :return:
        """
        conf = mol.GetConformer()
        d = np.array([])
        for fi, ci in mapping.items():
            fatom = self.followup.GetAtomWithIdx(fi)
            for neigh in fatom.GetNeighbors():
                ni = neigh.GetIdx()
                if ni not in mapping:
                    continue
                nci = mapping[ni]
                a = np.array(conf.GetAtomPosition(ci))
                b = np.array(conf.GetAtomPosition(nci))
                d = np.append(d, np.linalg.norm(a - b))
        return d

    def offness(self, mol: Chem.Mol, mapping: Dict[int, int]) -> float:
        """
        How many bonds are too long?
        :param mol:
        :param mapping:
        :return:
        """
        d = self.measure_map(mol, mapping)
        #return np.linalg.norm(d - 1.5)/(d.size*0.5) # 1.5 ang
        return sum(d > 2.5 ) * 3

