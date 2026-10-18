import csv
from collections import Counter

regions = {
    "N-terminal tail": (1, 30),
    "Fibril core": (31, 100),
    "C-terminal tail": (101, 140),
}

counts = Counter()

with open("results/residue_contacts.csv", newline="") as f:
    for row in csv.DictReader(f):
        r1 = int(row["Residue_1"])
        r2 = int(row["Residue_2"])

        region1 = next(
            name for name, (lo, hi) in regions.items()
            if lo <= r1 <= hi
        )
        region2 = next(
            name for name, (lo, hi) in regions.items()
            if lo <= r2 <= hi
        )

        pair = tuple(sorted((region1, region2)))
        counts[pair] += 1

print("Regional inter-chain contact counts")
print("-" * 50)

for pair, count in sorted(counts.items()):
    print(f"{pair[0]} vs {pair[1]}: {count}")

print("\nTotal contacts:", sum(counts.values()))
