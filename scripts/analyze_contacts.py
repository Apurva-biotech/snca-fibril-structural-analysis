import numpy as np
from Bio.PDB import PDBParser
from itertools import combinations

structure = PDBParser(QUIET=True).get_structure(
    "SNCA", "SNCA_fibril_2N0A.pdb"
)
chains = list(structure[0])

print("SNCA fibril: inter-chain C-alpha contacts")
print("Contact cutoff: 8.0 Angstrom")
print("-" * 65)
print("Pair | Minimum distance (A) | Contacts (<8 A)")

for c1, c2 in combinations(chains, 2):
    a1 = [r["CA"] for r in c1 if "CA" in r]
    a2 = [r["CA"] for r in c2 if "CA" in r]

    xyz1 = np.array([a.coord for a in a1])
    xyz2 = np.array([a.coord for a in a2])

    distances = np.linalg.norm(
        xyz1[:, None, :] - xyz2[None, :, :], axis=2
    )

    count = int(np.sum(distances < 8.0))
    minimum = float(np.min(distances))

    print(
        f"{c1.id}-{c2.id} | "
        f"{minimum:19.2f} | {count:14d}"
    )
