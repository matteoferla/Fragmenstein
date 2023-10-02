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

import logging
from rdkit import Chem
import numpy as np
from typing import List, Union, Sequence, Tuple, Optional, Dict
from rdkit import Chem
from rdkit.Chem import AllChem
import io, time
from openff.toolkit.topology import Molecule as OFFMolecule  # nomenclature idea copied from FEGrow
from openff.toolkit.topology import Topology as OFFTopology
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
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

    def init_pyrosetta(self):
        """
        This might be called if bad configuration is done?
        """
        self.journal.critical('Fritz does not use PyRosetta... Why was this called?')
        self.molar_energy_unit = mmu.kilocalorie_per_mole

    def __init__(self,
                 positioned_mol: Chem.Mol,
                 pdb_block: str,
                 resn: str = 'LIG',
                 resi: str = 1,
                 chain: str = 'X',
                 restrained_atomnames: Sequence[str] = (),
                 restraint_k: float = 1000.0,
                 mobile_radius: float = 8.0,
                 ):
        tick: float = time.time()
        self.resn: str = resn.strip()
        self.resi: int = int(resi)
        self.chain: str = chain
        # Igor does not have this as parameterisation is external
        self.positioned_mol: Chem.Mol = AllChem.AddHs(positioned_mol, addCoords=True)
        self.correct_pdbinfo(mol=self.positioned_mol, resn=self.resn, resi=self.resi, chain=self.chain)
        # this is apo
        self.pdb_block: str = pdb_block
        self.holo: mma.Modeller = self.plonk(pdb_block,
                                             self.positioned_mol)  # this is Fritz's plonk —Not Victor unlike Igor
        tock: float = time.time()
        self.journal.debug(f'Holo structure made {tock - tick}')
        self.simulation = self.create_simulation(restrained_atomnames=restrained_atomnames,
                                                 restraint_k=restraint_k,
                                                 mobile_radius=mobile_radius)
        self.unbound_simulation = self.create_simulation(restrained_atomnames=restrained_atomnames,
                                                         restraint_k=0,
                                                         mobile_radius=mobile_radius)
        self.shift_ligand(self.unbound_simulation, mm.Vec3(1_000, 0, 0) )
        tyck: float = time.time()
        self.journal.debug(f'Simulation created {tyck - tock}')

    from collections import defaultdict

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
            mol = self.positioned_mol
        if isinstance(apo, str):  # PDB block
            apo: mma.PDBFile = self.pdbblock_to_PDB(apo)
        # make a copy:
        holo: mma.Modeller = mma.Modeller(apo.topology, apo.positions)
        lig = OFFMolecule.from_rdkit(mol)
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
        holo.add(lig_topo, lig_pos)  # noqa mmu.Quantity is okay
        return holo

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

    def create_simulation(self,
                          restraint_k: float = 1_000,
                          restrained_atomnames: Sequence = (),
                          mobile_radius: float = 8.0,
                          forcefileds: Sequence[str] = ('amber14-all.xml', 'implicit/gbn2.xml')
                          ) -> mma.Simulation:
        """
        Creates a simulation object with the ligand harmonically constrained and the distal parts of the protein frozen.
        """
        # set up forcefield
        molecule = OFFMolecule.from_rdkit(self.positioned_mol)
        smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
        forcefield = mma.ForceField(*forcefileds)
        forcefield.registerTemplateGenerator(smirnoff.generator)
        # sort holo object which may have been tampered with...
        if isinstance(self.holo, mma.PDBFile):
            holo = mma.Modeller(self.holo.topology, self.holo.positions)
        self.holo.addHydrogens(forcefield, pH=7.0)
        # set up system
        system: mm.System = forcefield.createSystem(self.holo.topology,
                                                    nonbondedMethod=mma.NoCutoff,
                                                    nonbondedCutoff=1 * mmu.nanometer,
                                                    constraints=mma.HBonds)
        # restrain (harmonic constrain) the ligand
        if restraint_k:
            self.restrain(system, self.holo, k=restraint_k, atomnames=restrained_atomnames)
        integrator = mm.LangevinMiddleIntegrator(300 * mmu.kelvin, 1 / mmu.picosecond, 0.004 * mmu.picoseconds)
        simulation = mma.Simulation(self.holo.topology, system, integrator)
        simulation.context.setPositions(self.holo.positions)
        # freeze the distant parts of the protein
        self.freeze_distal(simulation, lig_resn=self.resn, radius=mobile_radius)
        return simulation

    def restrain(self, system: mm.System, pdb: Union[mma.PDBFile, mma.Modeller],
                 k: float = 1_000.0,
                 atomnames: Sequence[str] = (),
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
        hydrogen = mma.element.Element.getBySymbol('H')
        atomnames: List[str] = [n.strip() for n in atomnames]
        for atom in topology.atoms():
            if atom.element == hydrogen:
                continue
            if atomnames and atom.name not in atomnames:
                continue
            if atom.residue.name == 'LIG':
                restraint.addParticle(atom.index, positions[atom.index])
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
        # ## Freeze
        # freezing by setting mass to zero
        n = 0
        for atom in topology.atoms():
            if atom.residue.index not in ex_neighbors:
                simulation.system.setParticleMass(atom.index, 0. * mmu.amu)  # Dalton
            else:
                n += 1
        self.journal.info(f'{n} atoms were not frozen')
        return ex_neighbors

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

    def reanimate(self) -> Dict[str, Union[mmu.Quantity, str]]:
        tick: float = time.time()
        self.simulation.minimizeEnergy()
        tock: float = time.time()
        self.journal.debug(f'Reanimation! Minimisation of bound in {tock-tick}s')
        self.unbound_simulation.minimizeEnergy()
        tyck: float = time.time()
        self.journal.debug(f'Reanimation! Minimisation of unbound in {tyck-tock}s')
        data: Dict[str, Union[mmu.Quantity, str]] = {
            **{f'{k}_bound': v for k, v in self.get_potentials(self.simulation).items()},
            **{f'{k}_unbound': v for k, v in self.get_potentials(self.unbound_simulation).items()},
            'minimized_pdb': self.to_pdbblock()  # self.simulation by default
        }
        data['binding_dG'] = data['total_bound'] - data['total_unbound'] - data['CustomExternalForce_bound']
        return data

    def get_force_by_name(self,
                              name: str,
                              simulation: Optional[mma.Simulation]=None) -> mm.Force:
        if simulation is None:
            simulation = self.simulation
        for f in simulation.system.getForces():
            if f.getName() == name:
                return f
        raise ValueError(f'No force named {name}')

    def get_potentials(self, simulation: Optional[mma.Simulation]=None):
        if simulation is None:
            simulation = self.simulation
        potentials = {}
        potentials['total'] = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        for force in simulation.system.getForces():
            force.setForceGroup(30)
            potentials[force.getName()] = simulation.context\
                                                    .getState(getEnergy=True, groups={30})\
                                                    .getPotentialEnergy()
            force.setForceGroup(0)
        return potentials

    def remove_potential_by_name(self, name: str, simulation: Optional[mma.Simulation]=None):
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

    def to_pdbhandle(self, filehandle: io.TextIOWrapper, simulation: Optional[mma.Simulation]=None):
        if simulation is None:
            simulation = self.simulation
        positions: mmu.Quantity = simulation.context.getState(getPositions=True).getPositions()
        mma.PDBFile.writeFile(simulation.topology, positions, filehandle)

    def to_pdbfile(self, filename: str = 'fragmenstein.pdb', simulation: Optional[mma.Simulation]=None):
        with open(filename, 'w') as fh:
            self.to_pdbhandle(fh, simulation)

    def to_pdbblock(self, simulation: Optional[mma.Simulation]=None):
        iostr = io.StringIO()
        self.to_pdbhandle(iostr, simulation)
        iostr.seek(0)
        return iostr.read()

    def to_mol(self):
        topo, positions = self.get_topo_pos(self.simulation)
        lig_vec3: List[mm.Vec3] = [positions[atom.index].value_in_unit(mmu.angstrom) for atom in
                                   self.simulation.topology.atoms() if atom.residue.name == self.resn]
        n: int = self.positioned_mol.GetNumAtoms()
        if n != len(lig_vec3):
            self.journal.critical('Number of atoms discrepancy!')
        conf = Chem.Conformer(n)
        for mol_i, v in enumerate(lig_vec3):
            conf.SetAtomPosition(mol_i, (v.x, v.y, v.z))
        mol = Chem.Mol(self.positioned_mol)
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        return mol
