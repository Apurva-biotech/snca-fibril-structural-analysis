import csv
from collections import Counter

core_start, core_end = 31, 100
position_counts = Counter()
amino_acids = {}

with open("results/residue_contacts.csv", newline="") as f:
    for row in csv.DictReader(f):
        r1 = int(row["Residue_1"])
        r2 = int(row["Residue_2"])

        if core_start <= r1 <= core_end:
            position_counts[r1] += 1
            amino_acids[r1] = row["AA_1"]

        if core_start <= r2 <= core_end:
            position_counts[r2] += 1
            amino_acids[r2] = row["AA_2"]

print("Top 15 sequence positions by contact participation")
print("-" * 55)
print(f"{'Position':<12}{'AA':<8}{'Count':<10}")

for position, count in position_counts.most_common(15):
    print(f"{position:<12}{amino_acids[position]:<8}{count:<10}")

with open("results/sequence_hotspots.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Residue_Position", "Amino_Acid", "Contact_Participation_Count"])
    for position, count in sorted(position_counts.items()):
        writer.writerow([position, amino_acids[position], count])

print("\nSaved: results/sequence_hotspots.csv")
