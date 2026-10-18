import csv
from collections import Counter

core_start, core_end = 31, 100
cutoff = 8.0
contacts = Counter()

with open("results/residue_contacts.csv", newline="") as f:
    for row in csv.DictReader(f):
        r1 = int(row["Residue_1"])
        r2 = int(row["Residue_2"])
        d = float(row["CA_distance_A"])

        if d >= cutoff:
            continue

        if core_start <= r1 <= core_end:
            contacts[(row["Chain_1"], r1, row["AA_1"])] += 1

        if core_start <= r2 <= core_end:
            contacts[(row["Chain_2"], r2, row["AA_2"])] += 1

print("Top 20 fibril-core residue contact counts")
print("-" * 60)
print(f"{'Chain':<8}{'Residue':<10}{'AA':<8}{'Contacts':<10}")

for (chain, residue, aa), count in contacts.most_common(20):
    print(f"{chain:<8}{residue:<10}{aa:<8}{count:<10}")

print("\nUnique chain-specific core residues:", len(contacts))
