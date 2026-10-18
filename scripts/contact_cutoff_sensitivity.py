from Bio.PDB import PDBParser
from itertools import combinations
from math import dist

pdb_file = "SNCA_fibril_2N0A.pdb"
cutoffs = (6.0, 8.0, 10.0)

structure = PDBParser(QUIET=True).get_structure("SNCA", pdb_file)
chains = list(structure.get_chains())
counts = {cutoff: 0 for cutoff in cutoffs}

for chain1, chain2 in combinations(chains, 2):
    ca1 = [r["CA"].coord for r in chain1 if r.has_id("CA")]
    ca2 = [r["CA"].coord for r in chain2 if r.has_id("CA")]

    for xyz1 in ca1:
        for xyz2 in ca2:
            distance = dist(xyz1, xyz2)
            for cutoff in cutoffs:
                if distance < cutoff:
                    counts[cutoff] += 1

print("Inter-chain C-alpha contact sensitivity")
print("-" * 45)
for cutoff in cutoffs:
    print(f"Below {cutoff:.1f} Å: {counts[cutoff]} contacts")
