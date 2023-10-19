########################################################################################################################

__doc__ = \
    """
    See GitHub documentation
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2022 A.D."
__license__ = "MIT"
__citation__ = ""

import functools
import logging
from rdkit import Chem
import numpy as np
from typing import List, Union, Sequence, Tuple, Optional, Dict, Annotated
from rdkit import Chem
from rdkit.Chem import AllChem
import io, time
# nomenclature idea copied from FEGrow:
from openff.toolkit.topology import Molecule as OFFMolecule  # noqa
from openff.toolkit.topology import Topology as OFFTopology  # noqa
from openmmforcefields.generators import SMIRNOFFTemplateGenerator  # noqa
import openmm as mm
import openmm.app as mma
import openmm.unit as mmu
from pathlib import Path
from collections import defaultdict
import time


class Fritz:
    """
    Fritz is a helper class for Victor for use with OpenMM.
    It replaces Igor, the pyrosetta one.

    The two assistants are utterly different and share no attributes.
    """
    journal = logging.getLogger('Fragmenstein')
    forcefield_names = ('amber14-all.xml', 'implicit/gbn2.xml')
    molar_energy_unit = mmu.kilocalorie_per_mole
    integrator = mm.LangevinMiddleIntegrator(300 * mmu.kelvin, 1 / mmu.picosecond, 0.004 * mmu.picoseconds)

    def init_pyrosetta(self):
        """
        This might be called if bad configuration is done?
        """
        self.journal.critical('Fritz does not use PyRosetta... Why was this called?')

    def __init__(self,
                 prepped_mol: Chem.Mol,
                 pdb_block: str,
                 resn: str = 'LIG',
                 resi: str = 1,
                 chain: str = 'X',
                 restraining_atom_indices: Sequence[int] = (),
                 restraint_k: float = 1000.0,
                 mobile_radius: float = 8.0,
                 ):
        tick: float = time.time()
        self.resn: str = resn.strip()
        self.resi: int = int(resi)
        self.chain: str = chain
        # self.prepped_mol is the "initial" mol
        # Igor does not have this as parameterisation is external
        # self.prepped_mol is returned from `preminimized_undummied_mol`
        # so is technically not Victor.Monster.positioned_mol
        self.prepped_mol: Chem.Mol = AllChem.AddHs(prepped_mol, addCoords=True)
        self.correct_pdbinfo(mol=self.prepped_mol, resn=self.resn, resi=self.resi, chain=self.chain)
        # this is apo
        self.pdb_block: str = pdb_block
        self.holo: mma.Modeller = self.plonk(pdb_block,
                                             self.prepped_mol)  # this is Fritz's plonk —Not Victor unlike Igor
        tock: float = time.time()
        self.journal.debug(f'Holo structure made {tock - tick}')
        self.neighboring_res_idxs: List[int] = []  # computed later
        self.forcefield: mma.ForceField = self.create_forcefield()
        self.simulation = self.create_simulation(model=self.holo,
                                                 restraining_atom_indices=restraining_atom_indices,
                                                 restraint_k=restraint_k,
                                                 mobile_radius=mobile_radius)
        self.apo_simulation = self.create_simulation(model=self.pdbblock_to_PDB(self.pdb_block),
                                                     restraining_atom_indices=[],
                                                     restraint_k=0)
        self.ideal_ligand_simulation = self.create_simulation(model=self.rdkit_to_openMM(self.ideal_mol),
                                                              restraining_atom_indices=[],
                                                              restraint_k=0,
                                                              frozen=False)
        # self.shift_ligand(self.unbound_simulation, mm.Vec3(1_000, 0, 0))
        tyck: float = time.time()
        self.journal.debug(f'Simulation created {tyck - tock}')

    @staticmethod
    def correct_pdbinfo(mol: Chem.Mol, resn: str, resi: int, chain: str):
        """
        In Igor, RDKit to Params fixes it.
        In RDKit MolToPDBBlock https://github.com/rdkit/rdkit/blob/32655f5365e363ce13bd6b28e2e9e2544f8680bd/Code/GraphMol/FileParsers/PDBWriter.cpp#L66
        it is hidden.
        """
        missing = []
        dejavu = []
        counter = defaultdict(int)
        # first pass
        for atom in mol.GetAtoms():
            info: Union[None, Chem.AtomPDBResidueInfo] = atom.GetPDBResidueInfo()
            if info is None:
                missing.append(atom)
                continue
            name: str = info.GetName()
            if name.strip() in dejavu:
                # invalid.
                missing.append(atom)
            else:
                dejavu.append(name.strip())
                atom.SetProp('molFileAlias', name)
                info.SetResidueName(resn)
                info.SetResidueNumber(resi)

        def namegen(symbol):
            counter[symbol] += 1
            name = f'{symbol: >2}{counter[symbol]: <2}'
            if name.strip() in dejavu:
                return namegen(symbol)
            dejavu.append(name.strip())
            return name

        for atom in missing:
            symbol = atom.GetSymbol()
            name = namegen(symbol)
            info = Chem.AtomPDBResidueInfo(atomName=name,
                                           residueName=str(resn),
                                           residueNumber=int(resi),
                                           chainId=str(chain))
            atom.SetPDBResidueInfo(info)
            atom.SetProp('molFileAlias', name)

    def plonk(self,
              apo: Union[str, mma.PDBFile, mma.Modeller],
              mol: Optional[Chem.Mol] = None) -> mma.Modeller:
        """
        Plonk the ligand into the apo structure.
        """
        if mol is None:
            mol = self.prepped_mol
        if isinstance(apo, str):  # PDB block
            apo: mma.PDBFile = self.pdbblock_to_PDB(apo)
        # make a copy:
        holo: mma.Modeller = mma.Modeller(apo.topology, apo.positions)
        lig: mma.Modeller = self.rdkit_to_openMM(mol)
        holo.add(lig.topology, lig.positions)  # noqa mmu.Quantity is okay
        return holo

    def rdkit_to_openMM(self, mol: Chem.Mol) -> mma.Modeller:
        # rdkit AssignStereochemistryFrom3D previously applied
        lig = OFFMolecule.from_rdkit(mol, allow_undefined_stereo=True)
        # minor corrections:
        # there is no need to fix via
        # lig_topo._chains[0]._residues[0].name = 'LIG'
        for r_a, m_a in zip(mol.GetAtoms(), lig.atoms):
            assert r_a.GetSymbol() == m_a.symbol, 'Discrepancy'
            r_name = r_a.GetPDBResidueInfo().GetName()
            m_a.name = r_name
        # convert and merge
        lig_topo: mma.Topology = OFFTopology.from_molecules([lig]).to_openmm()
        lig_pos: mmu.Quantity = lig.conformers[0].to_openmm()
        return mma.Modeller(lig_topo, lig_pos)  # noqa mmu.Quantity is okay

    @staticmethod
    def pdbblock_to_PDB(pdb_block: str) -> mma.PDBFile:
        """
        Read a PDB block (string) into a PDBFile
        """
        assert isinstance(pdb_block, str), 'pdb_block must be a string'
        assert len(pdb_block) > 0, 'pdb_block must be a non-empty string'
        iostr = io.StringIO(pdb_block)
        iostr.seek(0)
        pdb = mma.PDBFile(iostr)
        return pdb

    def create_forcefield(self, forcefield_names: Optional[Sequence[str]] = None) -> mma.ForceField:
        """set up forcefield"""
        molecule = OFFMolecule.from_rdkit(self.ideal_mol, allow_undefined_stereo=True)
        smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
        if forcefield_names is None:
            forcefield_names = self.forcefield_names
        forcefield = mma.ForceField(*forcefield_names)
        forcefield.registerTemplateGenerator(smirnoff.generator)
        return forcefield

    @functools.cached_property
    def ideal_mol(self):
        ideal = AllChem.AddHs(self.prepped_mol)
        ideal.RemoveAllConformers()
        AllChem.EmbedMolecule(ideal, enforceChirality=True)
        return ideal

    def create_simulation(self,
                          model: Union[mma.PDBFile, mma.Modeller, None] = None,
                          restraint_k: float = 1_000,
                          restraining_atom_indices: Sequence[int] = (),
                          mobile_radius: float = 8.0,
                          frozen: bool = True,
                          ) -> mma.Simulation:
        """
        Creates a simulation object with the ligand harmonically constrained and the distal parts of the protein frozen.
        """
        # deal with model
        if isinstance(model, mma.Modeller):
            pass
        elif isinstance(model, mma.PDBFile):
            model = mma.Modeller(model.topology, model.positions)
        elif isinstance(self.holo, mma.PDBFile):
            model = mma.Modeller(self.holo.topology, self.holo.positions)
        else:
            model = self.holo
        # add Hydrogens
        model.addHydrogens(self.forcefield, pH=7.0)
        # set up system
        system: mm.System = self.forcefield.createSystem(model.topology,
                                                         nonbondedMethod=mma.NoCutoff,
                                                         nonbondedCutoff=1 * mmu.nanometer,
                                                         constraints=mma.HBonds)
        # restrain (harmonic constrain) the ligand
        if restraint_k:
            self.restrain(system, model, k=restraint_k, atom_indices=restraining_atom_indices)
        simulation = mma.Simulation(model.topology, system, self.integrator)
        simulation.context.setPositions(model.positions)
        # freeze the distant parts of the protein
        if frozen and len(self.neighboring_res_idxs) == 0:
            self.journal.debug('Freezing...')
            self.neighboring_res_idxs: List[int] = self.freeze_distal(simulation,
                                                                      lig_resn=self.resn,
                                                                      radius=mobile_radius)
            self.journal.debug(f'Frozen all bar {self.neighboring_res_idxs}')
        elif frozen:
            self.journal.debug('Freezing prespecified')
            self.freeze(simulation=simulation,
                        unfrozen_resn_idxs=self.neighboring_res_idxs)
        else:
            self.journal.debug('No freezing')
        return simulation

    def restrain(self, system: mm.System, pdb: Union[mma.PDBFile, mma.Modeller],
                 k: float = 1_000.0,
                 atom_indices: Sequence[int] = (),
                 **args):
        """
        This needs to be set before the simulation is created. I dont know why.
        """
        positions: mmu.Quantity
        topology: mma.topology.Topology
        topology, positions = self.get_topo_pos(pdb)
        # https://github.com/openmm/openmm-cookbook/blob/main/notebooks/cookbook/Restraining%20Atom%20Positions.ipynb
        restraint = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
        # restraint = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        # simulation.system.addForce(restraint)
        system.addForce(restraint)
        restraint.addGlobalParameter('k', k * mmu.kilojoules_per_mole / mmu.nanometer)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')

        # positions = simulation.context.getState(getPositions=True).getPositions()
        # ## get offset
        # The order seems to not be altered, so an offset works.
        hydrogen = mma.element.Element.getBySymbol('H')
        mm_atom: mma.topology.Atom
        rd_atom: Chem.Atom
        for mm_atom in topology.atoms():
            if mm_atom.residue.name == self.resn:
                offset = mm_atom.index
                break
        else:
            raise ValueError(f'{self.resn} missing')
        for mm_atom in topology.atoms():
            if mm_atom.residue.name != self.resn:
                continue
            elif mm_atom.element == hydrogen:
                continue
            # atom.index is C-style sequential index, atom.id is PDB "index".
            elif atom_indices and mm_atom.index - offset not in atom_indices:
                continue
            rd_atom: Chem.Atom = self.prepped_mol.GetAtomWithIdx(mm_atom.index - offset)
            if rd_atom.HasProp('_x'):
                # original atom positions
                atomic_xyz: mmu.Quantity = mm.Vec3(float(rd_atom.GetDoubleProp('_x')),
                                                   float(rd_atom.GetDoubleProp('_y')),
                                                   float(rd_atom.GetDoubleProp('_z'))) \
                                           * mmu.angstrom
            else:
                # Unlikely... but some hack may be at play. As these are added by `preminimized_undummied_mol`
                self.journal.debug('No _x property. Using position.')
                atomic_xyz: mmu.Quantity = positions[mm_atom.index]
            restraint.addParticle(mm_atom.index, atomic_xyz)
        self.journal.debug(f'N particles restrained: {restraint.getNumParticles()}')

    @staticmethod
    def get_topo_pos(obj) -> Tuple[mma.topology.Topology, mmu.Quantity]:
        """
        Get the topology and positions from a PDBFile, Modeller or Simulation object
        """
        if isinstance(obj, mma.PDBFile) or isinstance(obj, mma.Modeller):
            pdb = obj
            positions = pdb.positions  # noqa mma.PDBFile.positions is a property
            topology = pdb.topology
        elif isinstance(obj, mma.Simulation):
            simulation = obj
            positions: mmu.Quantity = simulation.context.getState(getPositions=True).getPositions()
            topology: mma.topology.Topology = simulation.topology
        else:
            raise ValueError('obj must be of type PDBFile, Modeller or Simulation')
        return topology, positions

    def freeze_distal(self,
                      simulation: mma.Simulation,
                      lig_resn: Optional[str] = None,
                      radius: float = 8.0,
                      **args
                      ) -> List[int]:
        """
        "Constrain" in the crystallographic solving sense.
        This function freezes _everything_ if applied before simulation is created.
        I cannot figure out why.
        """
        if lig_resn is None:
            lig_resn = self.resn
        # constrain...
        positions: mmu.Quantity
        topology: mma.topology.Topology
        topology, positions = self.get_topo_pos(simulation)
        # ## Get ligand centroid
        ligand_residues: List[mma.topology.Residue] = [res for res in topology.residues() if res.name == lig_resn]
        assert len(ligand_residues) == 1, f'There can only be one ligand. N here: {len(ligand_residues)}'
        ligand_residue: mma.topology.Residue = ligand_residues[0]
        ligand_center: mmu.Quantity = self.get_centroid(ligand_residue, positions)
        # ## Find neighbours
        neighbors: List[int] = list({atom.residue.index for atom in topology.atoms() if
                                     self.distance(positions[atom.index], ligand_center, mmu.angstrom) <= radius})
        self.journal.info(f'unfrozen: {len(neighbors)} residues')
        ex_neighbors: List[int] = neighbors + [ligand_residue.index]
        self.freeze(simulation=simulation, unfrozen_resn_idxs=ex_neighbors)
        return ex_neighbors

    def freeze(self, simulation: mma.Simulation, unfrozen_resn_idxs: List[int]):
        """
        Freeze the distal parts of the protein,
        ie. the residues that are not in `unfrozen_resn_idxs`
        This is called by `freeze_distal`.
        """
        positions: mmu.Quantity
        topology: mma.topology.Topology
        topology, positions = self.get_topo_pos(simulation)
        # freezing by setting mass to zero
        n = 0
        for atom in topology.atoms():
            if atom.residue.index not in unfrozen_resn_idxs:
                simulation.system.setParticleMass(atom.index, 0. * mmu.amu)  # Dalton
            else:
                n += 1
        self.journal.info(f'{n} atoms were not frozen')
        return unfrozen_resn_idxs

    @staticmethod
    def get_centroid(residue: mma.topology.Residue, positions: mmu.Quantity) -> mmu.Quantity:
        center: mmu.Quantity = mm.Vec3(0, 0, 0) * mmu.nanometer
        for atom in residue.atoms():
            center += positions[atom.index]
        return center / len(list(residue.atoms()))

    @staticmethod
    def distance(a: mmu.Quantity, b: mmu.Quantity, unit: mmu.Unit) -> np.float64:
        """
        The output is a float, not a Quantity, but multiply by unit to get the Quantity
        """
        d: mm.Vec3 = a.value_in_unit(unit) - b.value_in_unit(unit)
        return np.sqrt(d.x ** 2 + d.y ** 2 + d.z ** 2)

    def shift_ligand(self, simulation: mma.Simulation, amount: Union[mmu.Quantity, mm.Vec3, float]):
        if isinstance(amount, float):
            amount: mmu.Quantity = mm.Vec3(amount, 0, 0) * mmu.angstrom
        elif isinstance(amount, mm.Vec3):
            amount: mmu.Quantity = amount * mmu.angstrom
        else:
            pass
        positions: mmu.Quantity = simulation.context.getState(getPositions=True).getPositions()
        for atom in simulation.topology.atoms():
            if atom.residue.name != self.resn:
                continue
            positions[atom.index] = positions[atom.index] + amount
        simulation.context.setPositions(positions)

    def reanimate(self,
                  tolerance=10 * mmu.kilocalorie_per_mole / (mmu.nano * mmu.meter),
                  maxIterations: int = 0
                  ) -> Dict[str, Union[mmu.Quantity, str]]:
        if isinstance(tolerance, (float, int)):
            tolerance = float(tolerance) * mmu.kilocalorie_per_mole / (mmu.nano * mmu.meter)
        tick: float = time.time()
        self.simulation.minimizeEnergy(tolerance=tolerance, maxIterations=int(maxIterations))
        tock: float = time.time()
        self.journal.debug(f'Reanimation! Minimisation of bound in {tock - tick}s')
        self.apo_simulation.minimizeEnergy(tolerance=tolerance, maxIterations=int(maxIterations))
        self.ideal_ligand_simulation.minimizeEnergy(tolerance=tolerance, maxIterations=int(maxIterations))
        tyck: float = time.time()
        self.journal.debug(f'Reanimation! Minimisation of unbound in {tyck - tock}s')
        # sort data...
        data: Dict[str, Union[mmu.Quantity, str]] = {
            **{f'{k}_bound': v for k, v in self.get_potentials(self.simulation).items()},
            **{f'{k}_apo': v for k, v in self.get_potentials(self.apo_simulation).items()},
            **{f'{k}_ideal': v for k, v in self.get_potentials(self.ideal_ligand_simulation).items()},
            'minimized_pdb': self.to_pdbblock()  # self.simulation by default
        }
        for k in list(data.keys()):
            if '_bound' in k:
                term = k.replace('_bound', '')
                data[f'{term}_unbound'] = data.get(f'{term}_apo', 0) + data.get(f'{term}_ideal', 0)
        data['binding_dG'] = data['total_bound'] - data['total_unbound'] - data.get('CustomExternalForce_bound', 0.)
        return data

    def get_force_by_name(self,
                          name: str,
                          simulation: Optional[mma.Simulation] = None) -> mm.Force:
        if simulation is None:
            simulation = self.simulation
        for f in simulation.system.getForces():
            if f.getName() == name:
                return f
        raise ValueError(f'No force named {name}')

    def get_potentials(self, simulation: Optional[mma.Simulation] = None):
        if simulation is None:
            simulation = self.simulation
        potentials = {}
        potentials['total'] = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        for force in simulation.system.getForces():
            force.setForceGroup(30)
            potentials[force.getName()] = simulation.context \
                .getState(getEnergy=True, groups={30}) \
                .getPotentialEnergy()
            force.setForceGroup(0)
        return potentials

    def remove_potential_by_name(self, name: str, simulation: Optional[mma.Simulation] = None):
        """
        This exists to help do ``.remove_potential_by_name('CustomExternalForce')``
        """
        if simulation is None:
            simulation = self.simulation
        for fi, f in enumerate(simulation.system.getForces()):
            if f.getName() != name:
                continue
            simulation.system.removeForce(fi)
            return True
        return False

    def alter_restraint(self, new_k: float):
        self.journal.warning('`alter_restraint` is experimental and may not work.')
        f: mm.Force = self.get_force_by_name('CustomExternalForce')
        # get the index of the constant named _k_
        for i in range(f.getNumGlobalParameters()):
            if f.getGlobalParameterName(i) == 'k':
                break
        else:
            raise ValueError
        v = f.getGlobalParameterDefaultValue(i)
        f.setGlobalParameterDefaultValue(i, new_k)
        f.updateParametersInContext(self.simulation.context)

    def to_pdbhandle(self, filehandle: io.TextIOWrapper, simulation: Optional[mma.Simulation] = None):
        if simulation is None:
            simulation = self.simulation
        positions: mmu.Quantity = simulation.context.getState(getPositions=True).getPositions()
        mma.PDBFile.writeFile(simulation.topology, positions, filehandle)

    def to_pdbfile(self, filename: str = 'fragmenstein.pdb', simulation: Optional[mma.Simulation] = None):
        with open(filename, 'w') as fh:
            self.to_pdbhandle(fh, simulation)

    def to_pdbblock(self, simulation: Optional[mma.Simulation] = None):
        iostr = io.StringIO()
        self.to_pdbhandle(iostr, simulation)
        iostr.seek(0)
        return iostr.read()

    def to_mol(self, simulation: Optional[mma.Simulation] = None):
        if simulation is None:
            simulation: mma.Simulation = self.simulation
        topo, positions = self.get_topo_pos(simulation)
        lig_vec3: List[mm.Vec3] = [positions[atom.index].value_in_unit(mmu.angstrom) for atom in
                                   simulation.topology.atoms() if atom.residue.name == self.resn]
        n: int = self.prepped_mol.GetNumAtoms()
        if n != len(lig_vec3):
            self.journal.critical('Number of atoms discrepancy!')
        conf = Chem.Conformer(n)
        for mol_i, v in enumerate(lig_vec3):
            conf.SetAtomPosition(mol_i, (v.x, v.y, v.z))
        mol = Chem.Mol(self.prepped_mol)
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        return mol

    @classmethod
    def minimize_template(cls, apo_block: str,
                          forcefield_names: Optional[Sequence[str]] = None):
        pdb = Fritz.pdbblock_to_PDB(apo_block)
        apo = mma.Modeller(pdb.topology, pdb.positions)
        if forcefield_names is None:
            forcefield_names = cls.forcefield_names
        forcefield = mma.ForceField(*forcefield_names)
        apo.addHydrogens(forcefield, pH=7.0)
        system: mm.System = forcefield.createSystem(apo.topology,
                                                    nonbondedMethod=mma.NoCutoff,
                                                    nonbondedCutoff=1 * mmu.nanometer,
                                                    constraints=mma.HBonds)
        simulation = mma.Simulation(apo.topology, system, cls.integrator)
        simulation.context.setPositions(apo.positions)
        simulation.minimizeEnergy()
        positions: mmu.Quantity = simulation.context.getState(getPositions=True).getPositions()
        iostr = io.StringIO()
        mma.PDBFile.writeFile(simulation.topology, positions, iostr)
        iostr.seek(0)
        return iostr.read()

    def get_unfrozen_pymol_selection(self):
        """
        This is a pymol selection string.
        It acts by residue index only as Chain is being weird:
        `r.chain.id` is `''`, but ought to be `'A'`.
        (Index is 0-based, id is PDB-based w/ gaps and all)
        """
        t, p = self.get_topo_pos(self.simulation)
        _sele = []
        for r in t.residues():
            if r.index in self.neighboring_res_idxs:
                _sele.append(str(r.id))
        return 'resi ' + '+'.join(_sele)
