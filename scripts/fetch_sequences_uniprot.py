"""A script to fetch the desired protein from UniProt database using REST API"""
import requests
import sys
import urllib.parse

if len(sys.argv) != 4:
    sys.exit(
        "ERROR: Expected 3 arguments:\n"
        "fetch_sequences_uniprot.py <gene> <taxon> <output_file>")

gene, taxon, output_file = sys.argv[1:]

# URL-encode gene and taxon to avoid invalid characters
gene_safe = urllib.parse.quote(gene)
taxon_safe = urllib.parse.quote(taxon)

# Construct UniProt REST API URL
url = f"https://rest.uniprot.org/uniprotkb/stream?query=gene_exact:{gene_safe}+AND+organism:{taxon_safe}&format=fasta"

try:
    response = requests.get(url, timeout=30)
    response.raise_for_status()
except requests.exceptions.RequestException as e:
    sys.exit(f"ERROR: UniProt request failed: {e}")

# If no sequences found, create empty file and warn
if not response.text.strip():
    print(f"WARNING: No UniProt sequences found for gene='{gene}' and taxon='{taxon}'.")
    with open(output_file, "w") as out_fasta:
        out_fasta.write("")  # create empty file
else:
    with open(output_file, "w") as out_fasta:
        out_fasta.write(response.text)
    print(f"Successfully retrieved UniProt sequences for gene='{gene}' and taxon='{taxon}'.")