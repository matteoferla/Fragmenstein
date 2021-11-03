from __future__ import annotations
# # original smaller function
# from collections import namedtuple
#
# Parts = namedtuple('Parts',
#                    ['step', 'headers', 'coordinates', 'connections', 'tails'],
#                    defaults=(0, [], [], [], [])
#                    )

# # The following does not work
#
# block = Chem.MolToPDBBlock(mol)
# mol = Chem.MolFromPDBBlock(block)
# for atom in mol.GetAtoms():
#     sn = atom.GetPDBResidueInfo().GetSerialNumber()
#     atom.GetPDBResidueInfo().SetSerialNumber(sn+10)
# block = Chem.MolToPDBBlock(mol)

from textwrap import wrap

class MinimalPDBParser:
    """
    This purpose build PDB parser simply fixes the serial numbers.
    The reason is that writing a custom 50 line class is easier that
    having biopython or other non-builtin requirement as a requirement
    Importing the PDB into RDKit is inadvisable.
    """

    def __init__(self, block: str, remove_water=False, remove_other_hetatms=False, ligname="LIG"):
        self.ligname = ligname
        self.remove_water = remove_water
        self.remove_other_hetatms = remove_other_hetatms
        self.to_preserve_heteroResname = ["LIG"]
        self.water_resnames = ["HOH", "OH", "H"]
        if not remove_water:
            self.to_preserve_heteroResname +=  self.water_resnames

        self.step = 0
        # step = 0 header unfinished, 1 coordinates finished, 2 connections finished.
        self.headers = []
        self.coordinates = []
        self.connections = []
        self.tails = []
        self.parse(block)

    def parse(self, block:str) -> None:
        # ---- parse -----------------------------------------
        def starts_with(xrow, name): return xrow.find(name) == 0

        for row in block.split('\n'):
            row = row.strip()
            if row == '':
                continue
            elif starts_with(row, 'ATOM') or starts_with(row, 'HETATM'):
                if starts_with(row, 'HETATM'):
                    resType = row[17:21].strip()
                    if self.remove_other_hetatms:
                        if resType not in self.to_preserve_heteroResname:
                            continue
                    if self.remove_water:
                        if resType in self.water_resnames:
                            continue
                self.step = 1
                self.coordinates.append(row)
            elif starts_with(row, 'CONECT'):
                self.step = 2
                self.connections.append(row)
            elif starts_with(row, 'TER') or starts_with(row, 'END'):
                continue
            elif self.step == 0:
                self.headers.append(row)
            elif starts_with(row, 'ANISOU'):
                continue
            elif self.step == 1:
                print(f'What is {row}?')
            elif self.step > 1:
                self.tails.append(row)
            else:
                raise SyntaxError('Impossible')

    def __str__(self):
        return '\n'.join(self.headers + self.coordinates + self.connections + ['END'] + self.tails)

    def get_serial(self, entry: str) -> int:
        # ATOM    588 11 - 14
        return int(entry[6:12].strip())

    def get_residue_index(self, entry: str) -> int:
        # https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
        # 23 - 26        Integer       resSeq       Residue sequence number.
        return int(entry[22:26].strip())

    def get_chain(self, entry: str) -> str:
        # 22             Character     chainID      Chain identifier.
        return entry[21].strip()

    def get_residue_name(self, entry: str) -> int:
        # https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
        # 18 - 20        Residue name  resName      Residue name.
        return entry[17:20].strip()

    def get_max_serial(self) -> int:
        # assuming ordered
        return self.get_serial(self.coordinates[-1])

    def set_serial(self, entry: str, value: int) -> None:
        new = f'{entry[:6]}{value: >5}{entry[11:]}'
        i = self.coordinates.index(entry)
        self.coordinates[i] = new

    def offset_serials(self, offset: int) -> None:
        for entry in self.coordinates:
            original_serial = self.get_serial(entry)
            self.set_serial(entry, original_serial + offset)

    def offset_connections(self, offset:int) -> None:
        for i, entry in enumerate(self.connections):
            self.connections[i] = 'CONECT' + ''.join([f'{int(x)+offset: >5}' for x in wrap(entry[7:], 5)])

    def append(self, other: MinimalPDBParser):
        """
        Add a second parser data to it. But only its coordinates and connections.
        """
        offset = self.get_max_serial()
        other.offset_serials(offset)
        other.offset_connections(offset)
        self.coordinates += other.coordinates
        self.connections += other.connections

    def has_residue_index(self, index:int, chain: str):
        for entry in self.coordinates:
            if self.get_residue_index(entry) == index and self.get_chain(entry) == chain:
                return True
        else:
            return False

    def has_residue_name(self, name: str):
        """
        residue name, resn 3-letters
        """
        for entry in self.coordinates:
            if self.get_residue_name(entry) == name:
                return True
        else:
            return False


