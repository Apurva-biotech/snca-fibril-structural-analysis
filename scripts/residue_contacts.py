import csv
import numpy as np
from Bio.PDB import PDBParser
from itertools import combinations

structure = PDBParser(QUIET=True).get_structure(
    "SNCA", "SNCA_fibril_2N0A.pdb"
)
chains = list(structure[0])

with open("results/residue_contacts.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Chain_1", "Residue_1", "AA_1",
        "Chain_2", "Residue_2", "AA_2", "CA_distance_A"
    ])

    for c1, c2 in combinations(chains, 2):
        for r1 in c1:
            if "CA" not in r1:
                continue
            for r2 in c2:
                if "CA" not in r2:
                    continue

                distance = r1["CA"] - r2["CA"]
                if distance < 8.0:
                    writer.writerow([
                        c1.id, r1.id[1], r1.resname,
                        c2.id, r2.id[1], r2.resname,
                        round(distance, 3)
                    ])

print("Saved results/residue_contacts.csv")
