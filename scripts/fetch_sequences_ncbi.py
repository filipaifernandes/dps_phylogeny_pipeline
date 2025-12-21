"""A script to fetch the desired protein from NCBI Protein database"""

from Bio import Entrez
import sys
import time

if len(sys.argv) != 6:
    sys.exit(
        "ERROR: Expected 5 arguments:\n"
        "fetch_sequences_ncbi.py <gene> <taxon> <email> <max_seqs> <output>")
gene, taxon, email, max_seqs, output = sys.argv[1:]

Entrez.email = email

search_handle = Entrez.esearch(
    db="protein",
    term=f"{gene}[Gene] AND {taxon}[Organism]",
    retmax=int(max_seqs))

search_results = Entrez.read(search_handle)
search_handle.close()

ids = search_results.get("IdList", [])

if not ids:
    sys.exit(
        f"ERROR: No sequences found for gene='{gene}' and taxon='{taxon}'."
    )


with open(output, "w") as out_fasta:
    for seq_id in ids:
        fetch_handle = Entrez.efetch(
            db="protein",
            id=seq_id,
            rettype="fasta",
            retmode="text"
        )
        out_fasta.write(fetch_handle.read())
        fetch_handle.close()

        #Respect NCBI rate limits
        time.sleep(0.34)

print(f"Successfully retrieved {len(ids)} sequences from NCBI.")