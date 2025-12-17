"""A script to fetch the desired protein from UniProt database using REST API"""
import requests
import sys

gene, taxon, output = sys.argv[1:]

url = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    f"query=gene:{gene}+AND+organism:{taxon}"
    "&format=fasta"
)

response = requests.get(url)

with open(output, "w") as out:
    out.write(response.text)