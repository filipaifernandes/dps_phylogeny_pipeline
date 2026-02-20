"""A script to fetch the desired protein from UniProt database using REST API"""
#!/usr/bin/env python3
#!/usr/bin/env python3
import requests
import sys
import urllib.parse

if len(sys.argv) != 4:
    sys.exit(
        "ERROR: Expected 3 arguments:\n"
        "fetch_sequences_uniprot.py <gene> <taxon> <output_file>")

gene, taxon, output_file = sys.argv[1:]

# Encode gene and taxon for URL
gene_safe = urllib.parse.quote(gene)
taxon_safe = urllib.parse.quote(taxon)  # NO quotes around taxon

# Construct URL with gene_exact and properly encoded taxon
url = (
    f"https://rest.uniprot.org/uniprotkb/stream?"
    f"query=gene_exact:{gene_safe}+AND+organism:{taxon_safe}"
    "&format=fasta"
)

try:
    response = requests.get(url, timeout=30)
    response.raise_for_status()
except requests.exceptions.RequestException as e:
    sys.exit(f"ERROR: UniProt request failed: {e}")

if not response.text.strip():
    print(f"WARNING: No UniProt sequences found for gene='{gene}' and taxon='{taxon}'.")
    with open(output_file, "w") as out_fasta:
        out_fasta.write("")  # create empty file
else:
    with open(output_file, "w") as out_fasta:
        out_fasta.write(response.text)
    print(f"Successfully retrieved UniProt sequences for gene='{gene}' and taxon='{taxon}'.")