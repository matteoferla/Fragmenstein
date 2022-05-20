from rdkit.Chem import rdmolops

import json
from typing import Optional, Dict, List, Tuple, Set, Unpack, Union  # noqa cf. .legacy monkeypatch
from warnings import warn
from collections import defaultdict
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Geometry.rdGeometry import Point3D

from ._make_chimera import _MonsterChimera


class _MonsterRefine(_MonsterChimera):

    def posthoc_refine(self,
                       scaffold: Chem.Mol,
                       indices: Optional[List[int]] = None) -> Chem.Mol:
        """
        Given a scaffold and a list of indices, refine the scaffold.

        :param scaffold:
        :param indices: if absent, use all atoms
        :return:
        """
        if indices is None:
            indices = list(range(scaffold.GetNumAtoms()))
        refined = Chem.RWMol(scaffold)
        refconf = refined.GetConformer()
        positions = defaultdict(list)  # coordinates
        equivalence = defaultdict(list)  # atom indices of hits.
        for h in self.hits:
            if h.GetProp('_Name') in self.unmatched:
                continue
            hc = h.GetConformer()
            for k, v in self.get_positional_mapping(scaffold, h).items():
                positions[k].append([hc.GetAtomPosition(v).x, hc.GetAtomPosition(v).y, hc.GetAtomPosition(v).z])
                equivalence[k].append(f'{h.GetProp("_Name")}.{v}')
        for i in range(scaffold.GetNumAtoms()):
            if i not in indices:
                continue
            elif len(positions[i]) == 0:
                refined.GetAtomWithIdx(i).SetDoubleProp('_Stdev', 0.)
                refined.GetAtomWithIdx(i).SetDoubleProp('_Max', 0.)
                refined.GetAtomWithIdx(i).SetProp('_Origin', 'none')
                # warn(f'Atom {i}  {scaffold.GetAtomWithIdx(i).GetSymbol}/{refined.GetAtomWithIdx(i).GetSymbol} '+ \
                #     'in scaffold that has no positions.')
            else:
                p = np.mean(np.array(positions[i]), axis=0).astype(float)
                # sd = np.mean(np.std(np.array(positions[i]), axis=0)).astype(float)
                ds = [np.linalg.norm(p - pi) for pi in positions[i]]
                sd = np.std(ds)
                md = np.max(ds)
                refined.GetAtomWithIdx(i).SetProp('_Origin', json.dumps(equivalence[i]))
                refined.GetAtomWithIdx(i).SetDoubleProp('_Stdev', sd)
                refined.GetAtomWithIdx(i).SetDoubleProp('_Max', md)
                if self.average_position:
                    refconf.SetAtomPosition(i, Point3D(p[0], p[1], p[2]))
        Chem.SanitizeMol(refined,
                         sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS +
                                     Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY,
                         catchErrors=True)
        return refined

    # currently used only by partial
    # but maybe useful for other things
    def propagate_alternatives(self, fewer: List[Chem.Mol]) -> int:
        """
        Given the alt atoms strored in the Chem.Atom property ``_AltSymbol`` try those
        """
        pt = Chem.GetPeriodicTable()
        new = 0
        for template in list(fewer):
            for i, atom in enumerate(template.GetAtoms()):
                if atom.HasProp('_AltSymbol'):
                    alt = Chem.Mol(template)
                    aa = alt.GetAtomWithIdx(i)
                    aa.SetAtomicNum(pt.GetAtomicNumber(atom.GetProp('_AltSymbol')))
                    aa.ClearProp('_AltSymbol')
                    atom.ClearProp('_AltSymbol')
                    fewer.append(alt)
                    new += 1
        return new

    def pretweak(self) -> None:
        """
        What if the fragments were prealigned slightly? Really bad things happen.
        Nothing currently uses this without user interverntion.

        :return:
        """
        warn('This method is unreliable. Do not use it')
        ref = self.hits[0]
        for target in self.hits[1:]:
            A2B = list(self.get_positional_mapping(target, ref, 0.5).items())
            if A2B:
                rdMolAlign.AlignMol(target, ref, atomMap=A2B, maxIters=500)
            else:
                warn(f'No overlap? {A2B}')

    def sample_new_conformation(self, random_seed=None):
        """This method is intended for Multivictor.
        It generates a new conformation based on different random seeds.
        """
        scaffold = self.modifications["chimera"]
        atom_map = self.get_atom_map_fromProp(scaffold)
        if random_seed is None:
            random_seed = self.random_seed
        new_mol = self.place_from_map(target_mol=Chem.Mol(self.initial_mol), template_mol=scaffold, atom_map=atom_map,
                                      random_seed=random_seed)

        merging_mode = getattr(self, "merging_mode", "off")
        if merging_mode == 'off':
            pass
        elif merging_mode == 'full':
            pass
        elif merging_mode == 'partial':
            self.partial_blending()  # noqa will be in the subclasse partial
        elif merging_mode == 'none_permissive' or merging_mode == 'permissive_none' or \
                merging_mode == 'none':
            new_mol = self.posthoc_refine(new_mol)
        else:
            valid_modes = ('full', 'partial', 'none', 'none_permissive', 'off')
            raise ValueError(
                f"Merging mode can only be {'| '.join(valid_modes)}, not '{merging_mode}'")

        self.positioned_mol = new_mol
        return new_mol
