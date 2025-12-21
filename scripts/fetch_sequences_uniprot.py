"""A script to fetch the desired protein from UniProt database using REST API"""
import requests
import sys


if len(sys.argv) != 4:
    sys.exit(
        "ERROR: Expected 3 arguments:\n"
        "fetch_sequences_uniprot.py <gene> <taxon> <output>")

gene, taxon, output = sys.argv[1:]

url = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    f"query=gene:{gene}+AND+organism:{taxon}"
    "&format=fasta")

try:
    response = requests.get(url, timeout=30)
    response.raise_for_status()
except requests.exceptions.RequestException as e:
    sys.exit(f"ERROR: UniProt request failed: {e}")

if not response.text.strip():
    sys.exit(
        f"ERROR: No UniProt sequences found for gene='{gene}' and taxon='{taxon}'."
    )
with open(output, "w") as out_fasta:
    out_fasta.write(response.text)

print("Successfully retrieved UniProt sequences.")