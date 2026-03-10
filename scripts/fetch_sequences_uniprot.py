"""A script to fetch the desired protein from UniProt database using REST API"""
import sys
import requests

if len(sys.argv) != 6:
    sys.exit("Usage: fetch_sequences_uniprot.py <gene> <organism> <email> <max_seqs> <output>")

gene, organism, email, max_seqs, outfile = sys.argv[1:]

url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_name:{organism}&format=fasta&size={max_seqs}"

r = requests.get(url)
if r.status_code != 200 or not r.text.strip():
    raise Exception(f"UniProt request failed for gene={gene}, organism={organism}")

with open(outfile, "w") as f:
    f.write(r.text)

print(f"Saved UniProt sequences to {outfile}")
