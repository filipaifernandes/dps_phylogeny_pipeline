import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO

# IMPORTANT for Docker
plt.switch_backend("Agg")

alignment_file = sys.argv[1]
output_matrix = sys.argv[2]
output_plot = sys.argv[3]

alignment = AlignIO.read(alignment_file, "fasta")
n = len(alignment)

names = [record.id for record in alignment]

def identity(seq1, seq2):
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches / len(seq1)

matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        matrix[i, j] = identity(alignment[i].seq, alignment[j].seq)

df = pd.DataFrame(matrix, index=names, columns=names)

df.to_csv(output_matrix)

plt.figure(figsize=(10, 8))
sns.heatmap(df, cmap="viridis")
plt.title("Dps Sequence Identity Heatmap")
plt.tight_layout()
plt.savefig(output_plot)
