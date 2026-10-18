from Bio.PDB import PDBParser
from itertools import combinations
from math import dist
import csv

pdb_file = "SNCA_fibril_2N0A.pdb"
cutoff = 8.0

structure = PDBParser(QUIET=True).get_structure("SNCA", pdb_file)
chains = list(structure.get_chains())
results = []

for chain1, chain2 in combinations(chains, 2):
    ca1 = [(r.id[1], r.resname, r["CA"].coord)
           for r in chain1 if r.has_id("CA")]
    ca2 = [(r.id[1], r.resname, r["CA"].coord)
           for r in chain2 if r.has_id("CA")]

    contacts = 0
    min_distance = float("inf")

    for _, _, xyz1 in ca1:
        for _, _, xyz2 in ca2:
            d = dist(xyz1, xyz2)
            min_distance = min(min_distance, d)
            if d < cutoff:
                contacts += 1

    results.append([
        chain1.id, chain2.id, contacts, round(min_distance, 3)
    ])

with open("results/interface_contacts.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Chain_1", "Chain_2",
        "Contacts_below_8A", "Minimum_CA_distance_A"
    ])
    writer.writerows(results)

print("Adjacent-chain interfaces")
print("-" * 45)
for c1, c2, count, mind in results:
    if ord(c2) == ord(c1) + 1:
        print(f"{c1}-{c2}: {count} contacts; minimum distance {mind} Å")

print("\nSaved: results/interface_contacts.csv")
