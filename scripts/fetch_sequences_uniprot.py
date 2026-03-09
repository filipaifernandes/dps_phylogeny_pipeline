"""A script to fetch the desired protein from UniProt database using REST API"""
import sys
import requests

gene = sys.argv[1]
organism = sys.argv[2]
outfile = sys.argv[3]

url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_name:{organism}&format=fasta"

r = requests.get(url)
if r.status_code != 200 or not r.text.strip():
    raise Exception("UniProt request failed")

with open(outfile, "w") as f:
    f.write(r.text)

print(f"Saved UniProt sequences to {outfile}")
