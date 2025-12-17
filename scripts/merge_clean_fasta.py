"""A script to merge FASTA files and clean them"""
from Bio import SeqIO
import sys

inputs = sys.argv[1:-1]
output = sys.argv[-1]

seen = set()

with open(output, "w") as out:
    for fasta in inputs:
        for record in SeqIO.parse(fasta, "fasta"):
            seq = str(record.seq)
            if seq not in seen:
                seen.add(seq)
                SeqIO.write(record, out, "fasta")