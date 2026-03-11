from Bio import SeqIO

for record in SeqIO.parse("input.fasta", "fasta"):
    trunc = record.seq[53:207]   # 54–207
    record.seq = trunc
    record.id = "Dps1_trunc"
    record.description = ""
    SeqIO.write(record, "output.fasta", "fasta")
