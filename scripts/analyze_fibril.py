import numpy as np
from Bio.PDB import PDBParser

structure = PDBParser(QUIET=True).get_structure(
    "SNCA", "SNCA_fibril_2N0A.pdb"
)

print("SNCA alpha-synuclein fibril analysis")
print("-" * 50)
print("Models:", len(list(structure)))
print("Chain | CA atoms | Mean CA position (x, y, z)")

for chain in structure[0]:
    atoms = [residue["CA"] for residue in chain if "CA" in residue]
    if atoms:
        coords = np.array([atom.coord for atom in atoms])
        center = coords.mean(axis=0)
        print(
            f"{chain.id:>5} | {len(atoms):>9} | "
            f"{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}"
        )
