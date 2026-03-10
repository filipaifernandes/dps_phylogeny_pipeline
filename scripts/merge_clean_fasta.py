"""A script to merge FASTA files and clean them"""
from Bio import SeqIO
import sys
import re

inputs = sys.argv[1:-1]
output = sys.argv[-1]

seen = set()

MIN_LEN = 150
MAX_LEN = 300

valid_aa = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")

with open(output, "w") as out:
    for fasta in inputs:
        for record in SeqIO.parse(fasta, "fasta"):

            seq = str(record.seq).upper()

            # Remove duplicates
            if seq in seen:
                continue

            # Remove sequences with ambiguous residues
            if not valid_aa.match(seq):
                continue

            # Length filtering
            if len(seq) < MIN_LEN or len(seq) > MAX_LEN:
                continue

            seen.add(seq)
            SeqIO.write(record, out, "fasta")
