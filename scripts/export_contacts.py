import csv
import numpy as np
from itertools import combinations
from Bio.PDB import PDBParser

structure = PDBParser(QUIET=True).get_structure(
    "SNCA", "SNCA_fibril_2N0A.pdb"
)
chains = list(structure[0])
n = len(chains)
contact_matrix = np.zeros((n, n), dtype=int)

for i, j in combinations(range(n), 2):
    a = [r["CA"].coord for r in chains[i] if "CA" in r]
    b = [r["CA"].coord for r in chains[j] if "CA" in r]

    distances = np.linalg.norm(
        np.array(a)[:, None, :] - np.array(b)[None, :, :],
        axis=2
    )
    count = int(np.sum(distances < 8.0))
    contact_matrix[i, j] = count
    contact_matrix[j, i] = count

labels = [c.id for c in chains]

with open("results/contact_matrix.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Chain"] + labels)
    for label, row in zip(labels, contact_matrix):
        writer.writerow([label] + row.tolist())

print("Saved results/contact_matrix.csv")

try:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(contact_matrix, interpolation="nearest")
    ax.set_xticks(range(n), labels=labels)
    ax.set_yticks(range(n), labels=labels)
    ax.set_xlabel("Chain")
    ax.set_ylabel("Chain")
    ax.set_title("SNCA fibril: inter-chain C-alpha contacts (<8 Å)")
    fig.colorbar(im, ax=ax, label="Residue pairs")
    fig.tight_layout()
    fig.savefig("results/contact_heatmap.png", dpi=300)
    plt.close(fig)
    print("Saved results/contact_heatmap.png")
except ImportError:
    print("Matplotlib is missing; CSV saved, but heatmap was not generated.")
