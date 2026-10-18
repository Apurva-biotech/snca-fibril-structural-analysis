import csv
import matplotlib.pyplot as plt

positions = []
counts = []

with open("results/sequence_hotspots.csv", newline="") as f:
    for row in csv.DictReader(f):
        position = int(row["Residue_Position"])
        if 31 <= position <= 100:
            positions.append(position)
            counts.append(int(row["Contact_Participation_Count"]))

plt.figure(figsize=(12, 5))
plt.plot(positions, counts, linewidth=1.5)
plt.xlabel("Alpha-synuclein residue position")
plt.ylabel("Contact-participation count")
plt.title("Inter-chain Contact Profile of the Alpha-Synuclein Fibril Core")
plt.xlim(31, 100)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("results/core_contact_profile.png", dpi=300)
plt.show()

print("Saved: results/core_contact_profile.png")
