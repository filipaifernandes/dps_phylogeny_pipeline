"""A script to fetch the desired protein from NCBI Protein database"""

from Bio import Entrez
import sys

gene, taxon, email, max_seqs, output = sys.argv[1:]

Entrez.email = email

#Search NCBI protein database
handle = Entrez.esearch(
    db="protein",
    term=f"{gene}[Gene] AND {taxon}[Organism]",
    retmax=int(max_seqs)
)

ids = Entrez.read(handle)["IdList"]

#Fetch sequences and write as FASTA file
with open(output, "w") as out:
    for seq_id in ids:
        fetch = Entrez.efetch(db="protein", id=seq_id, rettype="fasta")
        out.write(fetch.read())