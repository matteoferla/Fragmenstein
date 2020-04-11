########################################################################################################################

__doc__ = \
    """
    ...
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = ""

########################################################################################################################


from typing import Dict, Union, List, Optional
from warnings import warn

try:
    from IPython.display import SVG, display
except ImportError:
    warn('No Jupyter notebook installed. `.draw_nicely` will not work.')
    SVG = lambda *args, **kwargs: print('Install IPython...')
    display = lambda *args, **kwargs: print('Install IPython...')

try:
    import pymol2
except ImportError:
    warn('No Pymol module installed. `.make_pse` will not work.')
    pymol2 = None

import numpy as np
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign, rdmolops
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry.rdGeometry import Point3D



##################################################################

class Fragmenstein:
    """
    Given a RDKit molecule and a series of hits it makes a spatially stitched together version of the initial molecule based on the hits.
    The reason is to do place the followup compound to the hits as faithfully as possible regardless of the screaming forcefields.

    * ``.scaffold`` is the combined version of the hits (rdkit.Chem.Mol object).
    * ``.chimera`` is the combined version of the hits, but with differing atoms made to match the followup (rdkit.Chem.Mol object).
    * ``.positioned_mol`` is the desired output (rdkit.Chem.Mol object)

    Note, the hits have to be spatially aligned —i.e. extracted from crystal structures in bond form.

    ``.get_positional_mapping``, which works also as a class method, creates a dictionary of mol_A atom index to mol_B atom index
    based on distance (cutoff 2&Aring;) and not MCS.

    The code works in two broad steps, first a scaffold is made, which is the combination of the hits (by position).
    Then the followup is placed. It is not embedded with constraint embed, which requires the reference molecule to have a valid geometry.
    ``.scaffold`` and ``.chimera`` and ``.positioned_mol`` absolutely do not have this.
    Novel side chains are added by aligning an optimised conformer against the closest 3-4 reference atoms.
    Note that ``.initial_mol`` is not touched. ``.positioned_mol`` may have lost some custom properties, but the atom idices are the same.
    """

    def __init__(self, mol: Chem.Mol, hits: List[Chem.Mol], debug_draw: bool = False):
        # starting attributes
        self.initial_mol = mol #untouched.
        #Chem.RemoveHs(self.initial_mol)
        self.hits = hits  # list of hits
        self._debug_draw = debug_draw # Jupyter notebook only.
        self.unmatched = []
        # derived attributes
        self.scaffold = self.merge_hits()  # merger of hits
        self.chimera = self.make_chimera()  # merger of hits but with atoms made to match the to-be-aligned mol
        self.positioned_mol = self.place_followup()  # to-be-aligned is aligned!

    def merge_hits(self) -> Chem.Mol:
        """
        Recursively stick the hits together and average the positions.

        :return: the rdkit.Chem.Mol object that will fill ``.scaffold``
        """
        hits = sorted(self.hits, key=lambda h: h.GetNumAtoms(), reverse=True)
        scaffold = Chem.Mol(hits[0])
        save_for_later = []
        for fragmentanda in hits[1:]:
            try:
                pairs = self._fragment_pairs(scaffold, fragmentanda)
                for anchor_index in pairs:
                    scaffold = self.merge(scaffold, fragmentanda,
                                          anchor_index=anchor_index,
                                          attachment_details=pairs[anchor_index])
            except ConnectionError:
                save_for_later.append(fragmentanda)
        for fragmentanda in save_for_later:
            try:
                pairs = self._fragment_pairs(scaffold, fragmentanda)
                for anchor_index in pairs:
                    scaffold = self.merge(scaffold, fragmentanda,
                                          anchor_index=anchor_index,
                                          attachment_details=pairs[anchor_index])
            except ConnectionError:
                self.unmatched.append(fragmentanda.GetProp("_Name"))
                warn(f'Hit {fragmentanda.GetProp("_Name")} has no connections! Skipping!')
        refined = self.posthoc_refine(scaffold)
        return refined

    def merge(self, scaffold:Chem.Mol, fragmentanda:Chem.Mol, anchor_index:int, attachment_details:List[Dict]) -> Chem.Mol:
        for detail in attachment_details:
            attachment_index = detail['idx_F']  # fragmentanda attachment_index
            scaffold_attachment_index = detail['idx_S']
            bond_type = detail['type']
            f = Chem.FragmentOnBonds(fragmentanda,
                                     [fragmentanda.GetBondBetweenAtoms(anchor_index, attachment_index).GetIdx()],
                                     addDummies=False)
            frag_split = []
            fragmols = Chem.GetMolFrags(f, asMols=True, fragsMolAtomMapping=frag_split, sanitizeFrags=False)
            if self._debug_draw:
                print(frag_split)
            # Get the fragment of interest.
            ii = 0
            for mol_N, indices in enumerate(frag_split):
                if anchor_index in indices:
                    break
                ii += len(indices)
            else:
                raise Exception
            frag = fragmols[mol_N]
            frag_anchor_index = indices.index(anchor_index)
            if self._debug_draw:
                self.draw_nicely(frag)
            combo = Chem.RWMol(rdmolops.CombineMols(scaffold, frag))
            scaffold_anchor_index = frag_anchor_index + scaffold.GetNumAtoms()
            if self._debug_draw:
                print(scaffold_anchor_index, scaffold_attachment_index, anchor_index, scaffold.GetNumAtoms())
                self.draw_nicely(combo)
            combo.AddBond(scaffold_anchor_index, scaffold_attachment_index, bond_type)
            Chem.SanitizeMol(combo,
                             sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS +
                                         Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY,
                             catchErrors=True)
            if self._debug_draw:
                self.draw_nicely(combo)
            scaffold = combo
        return scaffold

    def _fragment_pairs(self, scaffold:Chem.Mol, fragmentanda:Chem.Mol) -> Dict[int, List[Dict]]:
        A2B_mapping = self.get_positional_mapping(scaffold, fragmentanda)
        get_key = lambda d, v: list(d.keys())[list(d.values()).index(v)]
        if len(A2B_mapping) == 0:
            raise ConnectionError
        uniques = set(range(fragmentanda.GetNumAtoms())) - set(A2B_mapping.values())
        categories = self._categorise(fragmentanda, uniques)
        pairs = categories['pairs']
        for p in pairs:  # pairs:Dict[List[Dict]]
            for pp in pairs[p]:
                pp['idx_F'] = pp['idx']  # less ambiguous: fragmentanda index
                pp['idx_S'] = get_key(A2B_mapping, pp['idx'])  # scaffold index
        return pairs

    @classmethod
    def get_positional_mapping(cls, mol_A: Chem.Mol, mol_B: Chem.Mol, cutoff=2) -> Dict[int, int]:
        """
        Returns a map to convert overlapping atom of A onto B
        Cutoff 2 &Aring;.

        :param mol_A: first molecule (Chem.Mol) will form keys
        :param mol_B: second molecule (Chem.Mol) will form values
        :return: dictionary mol A atom idx -> mol B atom idx.
        """
        mols = [mol_A, mol_B]
        confs = [mols[i].GetConformers()[0] for i in range(len(mols))]
        distance = lambda a, b: ((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2) ** 0.5
        m = []
        for i in range(mols[0].GetNumAtoms()):
            v = []
            for j in range(mols[1].GetNumAtoms()):
                d = distance(confs[0].GetAtomPosition(i), confs[1].GetAtomPosition(j))
                v.append(d)
            m.append(v)
        dm = np.array(m)
        ## find the closest
        mapping = {}
        while 1 == 1:
            d = np.amin(dm)
            if d > cutoff:
                break
            w = np.where(dm == d)
            f, s = w[0][0], w[1][0]
            mapping[int(f)] = int(s)  # np.int64 --> int
            dm[f, :] = np.ones(dm.shape[1]) * 999
            dm[:, s] = np.ones(dm.shape[0]) * 999
        return mapping

    def _categorise(self, mol:Chem.Mol, uniques:set) -> Dict[str, Union[set, Dict]]:
        """
        What do the novel atoms do in terms of connectivity.
        Complicated dict output. Really ought to be SetProp of the atoms.

        * ``uniques`` are set of atoms to classify on
        * ``internals`` are unique atoms that are connected solely to unique atoms
        * ``attachments`` are non-unique atoms to which a unique atom connects
        * ``pairs`` is a dict of unique atom idx --> dict of ``idx`` --> attachment idx and ``type`` bond type.

        :param mol: molecule to describe
        :param uniques: set of indices that are new to this molecule
        :return:
        """
        #
        pairs = {}
        internals = set()
        attachments = set()
        for i in uniques:
            unique_atom = mol.GetAtomWithIdx(i)
            neighbours = {n.GetIdx() for n in unique_atom.GetNeighbors()}
            if len(neighbours - uniques) == 0:
                internals.add(i)
            else:
                i_attached = neighbours - uniques
                attachments |= i_attached
                pairs[i] = [{'idx': j,
                             'type': mol.GetBondBetweenAtoms(i, j).GetBondType()} for j in i_attached]
        anchors = uniques - internals
        if self._debug_draw:
            high = list(internals) + list(attachments) + list(anchors)
            color = {**{i: (0, 0.8, 0) for i in internals},
                     **{i: (0, 0, 0.8) for i in attachments},
                     **{i: (0.8, 0, 0.8) for i in anchors}}
            self.draw_nicely(mol, highlightAtoms=high, highlightAtomColors=color)
        return dict(uniques=uniques,
                    internals=internals,
                    attachments=attachments,
                    pairs=pairs
                    )

    def posthoc_refine(self, scaffold):
        """
        Averages the overlapping atoms.

        :param scaffold:
        :return:
        """
        refined = Chem.RWMol(scaffold)
        refconf = refined.GetConformer()
        positions = defaultdict(list)
        for h in self.hits:
            hc = h.GetConformer()
            for k, v in self.get_positional_mapping(scaffold, h).items():
                positions[k].append([hc.GetAtomPosition(v).x, hc.GetAtomPosition(v).y, hc.GetAtomPosition(v).z])
        for i in range(scaffold.GetNumAtoms()):
            if len(positions[i]) == 0:
                continue
                warn(f'Atom {i}  {scaffold.GetAtomWithIdx(i).GetSymbol}/{refined.GetAtomWithIdx(i).GetSymbol} in scaffold that has no positions.')
            p = np.mean(np.array(positions[i]), axis=0).astype(float)
            refconf.SetAtomPosition(i, Point3D(p[0], p[1], p[2]))
        Chem.SanitizeMol(refined,
                         sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS +
                                     Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY,
                         catchErrors=True)
        return refined

    def make_chimera(self) -> Chem.Mol:
        """
        This is to avoid corner cases. E.g. here the MCS is ringMatchesRingOnly=True and AtomCompare.CompareAny,
        while for the positioning this is not the case.

        :return:
        """
        mcs = rdFMCS.FindMCS([self.scaffold, self.initial_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareAny,
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        scaffold_match = self.scaffold.GetSubstructMatch(common)
        followup_match = self.initial_mol.GetSubstructMatch(common)
        atomMap = [(followup_at, scaffold_at) for followup_at, scaffold_at in zip(followup_match, scaffold_match)]
        assert followup_match, 'No matching structure? All dummy atoms'
        if self._debug_draw:
            self.draw_nicely(common)
        ## make the scaffold more like the followup to avoid weird matches.
        chimera = Chem.RWMol(self.scaffold)
        for i in range(common.GetNumAtoms()):
            if common.GetAtomWithIdx(i).GetSymbol() == '*':  # dummies.
                v = {'F': 1, 'Br': 1, 'Cl': 1, 'B': 3, 'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 6}
                wanted = self.initial_mol.GetAtomWithIdx(followup_match[i])
                owned = self.scaffold.GetAtomWithIdx(scaffold_match[i])
                diff_valance = owned.GetExplicitValence() - v[wanted.GetSymbol()]
                if wanted.GetSymbol() in ('F', 'Br', 'Cl', 'C', 'H') and diff_valance > 0:
                    continue # cannot change this.
                elif owned.GetExplicitValence() > 4 and wanted.GetSymbol() not in ('P', ):
                    continue
                else:
                    chimera.ReplaceAtom(scaffold_match[i], Chem.Atom(wanted))
                    if diff_valance > 0:
                        chimera.GetAtomWithIdx(scaffold_match[i]).SetFormalCharge(diff_valance)
        try:
            chimera.UpdatePropertyCache()
        except Chem.AtomValenceException as err:
            warn('Valance issue'+str(err))
        return chimera

    @property
    def num_common(self) -> int:
        mcs = rdFMCS.FindMCS([self.scaffold, self.initial_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder)
        return Chem.MolFromSmarts(mcs.smartsString).GetNumAtoms()

    @property
    def percent_common(self) -> int:
        return round(self.num_common/self.initial_mol.GetNumAtoms()*100)

    def place_followup(self, mol:Chem.Mol=None) -> Chem.Mol:
        # Note none of this malarkey: AllChem.MMFFOptimizeMolecule(ref)
        # prealignment
        if mol is None:
            mol = self.initial_mol
        sextant = Chem.Mol(mol)
        Chem.SanitizeMol(sextant)
        AllChem.EmbedMolecule(sextant)
        AllChem.MMFFOptimizeMolecule(sextant)
        mcs = rdFMCS.FindMCS([self.chimera, mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder)
        common = Chem.MolFromSmarts(mcs.smartsString)
        chimera_match = self.chimera.GetSubstructMatch(common)
        followup_match = mol.GetSubstructMatch(common)
        atomMap = [(followup_at, chimera_at) for followup_at, chimera_at in zip(followup_match, chimera_match)]
        assert followup_match, 'No matching structure? All dummy atoms'
        rdMolAlign.AlignMol(sextant, self.chimera, atomMap=atomMap, maxIters=500)
        if self._debug_draw:
            self.draw_nicely(mol)
            self.draw_nicely(self.chimera)
            self.draw_nicely(common)
            print('followup/probe/mobile/candidate', followup_match)
            print('chimera/ref/scaffold', chimera_match)

        putty = Chem.Mol(sextant)
        pconf = putty.GetConformer()
        chimera_conf = self.chimera.GetConformer()
        pd = dict(atomMap)
        uniques = set()
        for i in range(putty.GetNumAtoms()):
            if i in pd:
                pconf.SetAtomPosition(i, chimera_conf.GetAtomPosition(pd[i]))
            else:
                uniques.add(i)
        # we be using a sextant for dead reckoning!
        categories = self._categorise(sextant, uniques)
        if self._debug_draw:
            print('internal', categories['internals'])
        for unique_idx in categories['pairs']:  # attachment unique indices
            sights = set()
            for pd in categories['pairs'][unique_idx]:
                first_sight = pd['idx']
                sights.add((first_sight, first_sight))
                neighs = [i.GetIdx() for i in sextant.GetAtomWithIdx(first_sight).GetNeighbors() if
                          i.GetIdx() not in uniques]
                for n in neighs:
                    sights.add((n, n))
            team = self._recruit_team(mol, unique_idx, categories)
            rdMolAlign.AlignMol(sextant, putty, atomMap=list(sights), maxIters=500)
            sconf = sextant.GetConformer()
            if self._debug_draw:
                print(f'alignment atoms for {unique_idx} ({team}): {sights}')
            for atom_idx in team:
                pconf.SetAtomPosition(atom_idx, sconf.GetAtomPosition(atom_idx))

        AllChem.SanitizeMol(putty)
        return putty

    def _recruit_team(self, mol:Chem.Mol, starting:set, categories:dict, team:Optional[set]=None) -> set:
        if team is None:
            team = set()
        team.add(starting)
        for atom in mol.GetAtomWithIdx(starting).GetNeighbors():
            i = atom.GetIdx()
            if i in categories['internals'] and i not in team:
                team = self._recruit_team(mol, i, categories, team)
        return team

    def make_pse(self, filename='test.pse'):
        assert '.pse' in filename, 'Must be a pymol pse extension!'
        with pymol2.PyMOL() as pymol:
            tints = iter(['wheat', 'palegreen', 'lightblue', 'paleyellow', 'lightpink', 'palecyan', 'lightorange', 'bluewhite'])
            #pymol.cmd.bg_color('white')
            for h, hit in enumerate(self.hits):
                pymol.cmd.read_molstr(Chem.MolToMolBlock(hit, kekulize=False), f'hit{h}')
                pymol.cmd.color(next(tints), f'hit{h} and name C*')
            pymol.cmd.read_molstr(Chem.MolToMolBlock(self.scaffold, kekulize=False), f'scaffold')
            pymol.cmd.color('tv_blue', f'scaffold and name C*')
            pymol.cmd.read_molstr(Chem.MolToMolBlock(self.chimera, kekulize=False), f'chimera')
            pymol.cmd.color('cyan', f'chimera and name C*')
            pymol.cmd.read_molstr(Chem.MolToMolBlock(self.positioned_mol, kekulize=False), f'followup')
            pymol.cmd.color('tv_green', f'followup and name C*')
            pymol.cmd.hide('sticks')
            pymol.cmd.hide('cartoon') # there should not be....
            pymol.cmd.show('lines', 'not polymer')
            pymol.cmd.show('sticks', 'followup or chimera')
            pymol.cmd.save(filename)

    def draw_nicely(self, mol, show=True, **kwargs):
        """
        Draw with atom indices for Jupyter notebooks.


        :param mol:
        :param kwargs: Key value pairs get fed into ``PrepareAndDrawMolecule``.
        :return:
        """
        d = rdMolDraw2D.MolDraw2DSVG(400, 400)
        d.drawOptions().addAtomIndices = True
        d.drawOptions().addStereoAnnotation = True
        x = Chem.Mol(mol)
        AllChem.Compute2DCoords(x)
        rdMolDraw2D.PrepareAndDrawMolecule(d, x, **kwargs)
        d.FinishDrawing()
        if show:
            display(SVG(d.GetDrawingText()))
        return d

    def save_commonality(self, filename:Optional[str]=None):
        mcs = rdFMCS.FindMCS([self.chimera, self.positioned_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        match = self.positioned_mol.GetSubstructMatch(common)
        d = self.draw_nicely(self.positioned_mol, show=False, highlightAtoms=match)
        if filename is None:
            return d
        else:
            assert '.svg' in filename, 'Can only save SVGs.'
            with open(filename, 'w') as w:
                w.write(d.GetDrawingText())


    def pretweak(self) -> None:
        """
        What if the fragments were prealigned slightly? Really bad things.

        :return:
        """
        print('This method is unreliable. Do not use it')
        ref = self.hits[0]
        for target in self.hits[1:]:
            A2B = list(self.get_positional_mapping(target, ref, 0.5).items())
            if A2B:
                rdMolAlign.AlignMol(target, ref, atomMap=A2B, maxIters=500)
            else:
                print(f'No overlap? {A2B}')

def test():
    hits = [Chem.MolFromMolFile(f'../Mpro/Mpro-{i}_0/Mpro-{i}_0.mol') for i in ('x0692', 'x0305', 'x1249')]
    followup = Chem.MolFromSmiles('CCNc1nc(CCS)c(C#N)cc1CN1C(CCS)CN(C(C)=O)CC1')
    Fragmenstein(followup, hits).make_pse('test.pse')

if __name__ == '__main__':
    test()
