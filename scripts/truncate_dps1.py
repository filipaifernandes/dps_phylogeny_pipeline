import sys
from Bio import SeqIO

input_file = sys.argv[1]
output_file = sys.argv[2]

truncated_records = []

for i, record in enumerate(SeqIO.parse(input_file, "fasta"), start=1):
    trunc_seq = record.seq[53:207]  # 54–207
    record.seq = trunc_seq
    record.id = f"{record.id}_trunc"
    record.description = ""
    truncated_records.append(record)

SeqIO.write(truncated_records, output_file, "fasta")
