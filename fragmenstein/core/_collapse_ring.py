import json
from warnings import warn
from rdkit.Geometry.rdGeometry import Point3D
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional, Dict, List

class Ring:
    def __init__(self):
        pass

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

    def _print_stored(self, mol: Chem.Mol):
        for atom in mol.GetAtoms():
            try:
                if atom.GetIntProp('_ori_i') == -1:
                    print(atom.GetIdx(),
                          atom.GetSymbol(),
                          atom.GetIntProp('_ori_i'),
                          atom.GetProp('_ori_name'),
                          atom.GetProp('_ori_is'),
                          atom.GetProp('_neighbors')
                    )
                else:
                    print(atom.GetIdx(),
                          atom.GetSymbol(),
                          atom.GetIntProp('_ori_i'),
                          atom.GetProp('_ori_name'))
            except KeyError:
                print(atom.GetIdx(),
                      atom.GetSymbol(), '!!!!')

    def _get_ori_i(self, mol: Chem.Mol, include_collapsed=True):
        indices = [atom.GetIntProp('_ori_i') for atom in mol.GetAtoms()]
        if include_collapsed:
            for atom in self._get_collapsed_atoms(mol):
                indices.extend(json.loads(atom.GetProp('_ori_is')))
        else:
            pass
        return indices

    def _get_collapsed_atoms(self, mol: Chem.Mol) -> List[Chem.Atom]:
        return [atom for atom in mol.GetAtoms() if atom.GetIntProp('_ori_i') == -1]

    def _get_new_index(self, mol: Chem.Mol, old: int, search_collapsed=True, name_restriction:Optional[str]=None) -> int:
        """
        Given an old index check in ``_ori_i`` for what the current one is.
        NB. ring placeholder will be -1 and these also have ``_ori_is``. a JSON of the ori_i they summarise.

        :param mol:
        :param old: old index
        :param search_collapsed: seach also in ``_ori_is``
        :parm name_restriction: restrict to original name.
        :return:
        """
        for i, atom in enumerate(mol.GetAtoms()):
            if name_restriction is not None and atom.GetProp('_ori_name') != name_restriction:
                pass
            elif atom.GetIntProp('_ori_i') == old:
                return i
            elif search_collapsed and \
                    atom.HasProp('_ori_is') and \
                    old in json.loads(atom.GetProp('_ori_is')):
                return i
            else:
                pass
        else:
            raise ValueError(f'Cannot find {old}')

    def _get_mystery_ori_i(self, mol: Chem.Mol):
        """
        Collapsed rings may remember neighbors that do not exist any more...

        :param mol:
        :return:
        """
        present = self._get_ori_i(mol)
        absent = []
        for atom in self._get_collapsed_atoms(mol):
            neighss = json.loads(atom.GetProp('_neighbors'))
            absent.extend([n for neighs in neighss for n in neighs if n not in present])
        return absent

    _collapsed_ring_offset=0
    def _offset_collapsed_ring(self, mol:Chem.Mol):
        """
        This is to prevent clashes.
        The numbers of the ori indices stored in collapsed rings are offset by the class variable (_collapsed_ring_offset)
        multiples of 100. (autoincrements to avoid dramas)

        :param mol:
        :return:
        """
        self._collapsed_ring_offset += 100
        for atom in self._get_collapsed_atoms(mol):
            old = json.loads(atom.GetProp('_ori_is'))
            new = [i+self._collapsed_ring_offset for i in old]
            atom.SetProp('_ori_is', json.dumps(new))
            old2new = dict(zip(old, new))
            print('UPDATE', old2new)
            old_neighss = json.loads(atom.GetProp('_neighbors'))
            new_neighss = [[old2new[i] if i in old2new else i for i in old_neighs] for old_neighs in old_neighss]
            atom.SetProp('_neighbors', json.dumps(new_neighss))

    def _offset_origins(self, mol:Chem.Mol):
        """
        This is to prevent clashes.

        :param mol:
        :return:
        """
        self._collapsed_ring_offset += 100
        old2new = {}
        for atom in mol.GetAtoms():
            if atom.GetIntProp('_ori_i') != -1:
                o = atom.GetIntProp('_ori_i')
                n = o + self._collapsed_ring_offset
                atom.SetIntProp('_ori_i', n)
                old2new[o] = n
        for atom in self._get_collapsed_atoms(mol):
            old = json.loads(atom.GetProp('_ori_is'))
            new = [i + self._collapsed_ring_offset for i in old]
            atom.SetProp('_ori_is', json.dumps(new))
            old2new = {**old2new, **dict(zip(old, new))}
            print('UPDATE', old2new)
            old_neighss = json.loads(atom.GetProp('_neighbors'))
            new_neighss = [[old2new[i] if i in old2new else i for i in old_neighs] for old_neighs in old_neighss]
            atom.SetProp('_neighbors', json.dumps(new_neighss))

    def _renumber_original_indices(self, mol: Chem.Mol,
                                   mapping: Dict[int, int],
                                   name_restriction:Optional[str]=None):
        """

        :param mol:
        :param mapping: old index to new index dictionary
        :param name_restriction:
        :return:
        """
        for atom in mol.GetAtoms():
            if name_restriction is not None and atom.HasProp('_ori_name') and atom.GetProp('_ori_name') != name_restriction:
                continue
            i = atom.GetIntProp('_ori_i')
            if i == -1:
                # original indices
                ori = json.loads(atom.GetProp('_ori_is'))
                alt = [dd if dd not in mapping else mapping[dd] for dd in ori]
                atom.SetProp('_ori_is', json.dumps(alt))
                # original neighbors
                ori = json.loads(atom.GetProp('_neighbors'))
                alt = [[dd if dd not in mapping else mapping[dd] for dd in inner] for inner in ori]
                atom.SetProp('_neighbors', json.dumps(alt))
            elif i in mapping:
                atom.SetIntProp('_ori_i', mapping[i])
            else:
                pass


    def expand_ring(self, mol: Chem.Mol):
        """
        Undoes collapse ring

        :param mol:
        :return:
        """
        mod = Chem.RWMol(mol)
        conf = mod.GetConformer()
        to_be_waited_for = []
        for atom in self._get_collapsed_atoms(mol):
            elements = json.loads(atom.GetProp('_elements'))
            neighbors = json.loads(atom.GetProp('_neighbors'))
            ori_i = json.loads(atom.GetProp('_ori_is'))
            xs = json.loads(atom.GetProp('_xs'))
            ys = json.loads(atom.GetProp('_ys'))
            zs = json.loads(atom.GetProp('_zs'))
            bonds = json.loads(atom.GetProp('_bonds'))
            new_is = []
            for i in range(len(elements)):
                n = mod.AddAtom(Chem.Atom(elements[i]))
                new_is.append(n)
                natom = mod.GetAtomWithIdx(n)
                conf.SetAtomPosition(n, Point3D(*[axis[i] for axis in (xs, ys, zs)]))
                natom.SetIntProp('_ori_i', ori_i[i])
                natom.SetDoubleProp('_x', xs[i])
                natom.SetDoubleProp('_y', ys[i])
                natom.SetDoubleProp('_z', zs[i])
            for i in range(len(elements)):
                new_i = new_is[i]
                for old_neigh, bond in zip(neighbors[i], bonds[i]):
                    bt = getattr(Chem.BondType, bond)
                    try:
                        new_neigh = self._get_new_index(mod, old_neigh, search_collapsed=False)
                        if not mod.GetBondBetweenAtoms(new_i, new_neigh):
                            mod.AddBond(new_i, new_neigh, bt)
                    except ValueError:
                        to_be_waited_for.append((old_neigh, bt))
        for old_neigh, bt in to_be_waited_for:
            try:
                new_neigh = self._get_new_index(mod, old_neigh, name_restriction=atom.GetProp('_ori_name'))
                if not mod.GetBondBetweenAtoms(new_i, new_neigh):
                    mod.AddBond(new_i, new_neigh, bt)
            except KeyError as err:
                warn(str(err))
        for a in reversed(range(mod.GetNumAtoms())):
            if mod.GetAtomWithIdx(a).GetIntProp('_ori_i') == -1:
                mod.RemoveAtom(a)
        return mod.GetMol()

    def collapse_ring(self, mol: Chem.Mol):
        """
        Collapses a ring(s) into a single dummy atom(s).
        Stores data as JSON in the atom.

        :param mol:
        :return:
        """
        self.store_positions(mol)
        mod = Chem.RWMol(mol)
        conf = mod.GetConformer()
        center_idxs = []
        morituri = []
        old2center = defaultdict(list)
        for atomset in mol.GetRingInfo().AtomRings():
            morituri.extend(atomset)
            neighs = []
            neighbonds = []
            bonds = []
            xs = []
            ys = []
            zs = []
            elements = []
            # add elemental ring
            c = mod.AddAtom(Chem.Atom('C'))
            center_idxs.append(c)
            central = mod.GetAtomWithIdx(c)
            name = mol.GetProp('_Name') if mol.HasProp('_Name') else '???'
            central.SetProp('_ori_name', name),
            # get data for storage
            for i in atomset:
                old2center[i].append(c)
                atom = mol.GetAtomWithIdx(i)
                neigh_i = [a.GetIdx() for a in atom.GetNeighbors()]
                neighs.append(neigh_i)
                bond = [mol.GetBondBetweenAtoms(i, j).GetBondType().name for j in neigh_i]
                bonds.append(bond)
                pos = conf.GetAtomPosition(i)
                xs.append(pos.x)
                ys.append(pos.y)
                zs.append(pos.z)
                elements.append(atom.GetSymbol())
            # store data in elemental ring
            central.SetIntProp('_ori_i', -1)
            central.SetProp('_ori_is', json.dumps(atomset))
            central.SetProp('_neighbors', json.dumps(neighs))
            central.SetProp('_xs', json.dumps(xs))
            central.SetProp('_ys', json.dumps(ys))
            central.SetProp('_zs', json.dumps(zs))
            central.SetProp('_elements', json.dumps(elements))
            central.SetProp('_bonds', json.dumps(bonds))
            conf.SetAtomPosition(c, Point3D(*[sum(axis) / len(axis) for axis in (xs, ys, zs)]))
        for atomset, center_i in zip(mod.GetRingInfo().AtomRings(), center_idxs):
            # bond to elemental ring
            central = mod.GetAtomWithIdx(center_i)
            neighss = json.loads(central.GetProp('_neighbors'))
            bondss = json.loads(central.GetProp('_bonds'))
            for neighs, bonds in zip(neighss, bondss):
                for neigh, bond in zip(neighs, bonds):
                    if neigh not in atomset:
                        bt = getattr(Chem.BondType, bond)
                        if neigh not in morituri:
                            mod.AddBond(center_i, neigh, bt)
                        else:
                            for other_center_i in old2center[neigh]:
                                if center_i != other_center_i:
                                    if not mod.GetBondBetweenAtoms(center_i, other_center_i):
                                        mod.AddBond(center_i, other_center_i, bt)
                                    break
                            else:
                                raise ValueError(f'Cannot find what {neigh} became')
        for i in sorted(set(morituri), reverse=True):
            mod.RemoveAtom(self._get_new_index(mod, i))
        return mod.GetMol()


if __name__ == '__main__':
    pass
    # for m in (Chem.MolFromSmiles('c1ccccc1CCCCN'), #benzene
    #           Chem.MolFromSmiles('c1ccccc1CCCCNCc2ccccc2'), #two benzenes
    #           Chem.MolFromSmiles('c1ccc2ccccc2c1CC')): # napthalene.
    #     AllChem.Compute2DCoords(m)
    #     Ring().store_positions(m)
    #     print('original')
    #     display(m)
    #     mod = Ring().collapse_ring(m)
    #     print('collapsed')
    #     display(mod)
    #     print('expanded')
    #     mod = Ring().expand_ring(mod)
    #     display(mod)