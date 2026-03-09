"""A script to fetch the desired protein from UniProt database using REST API"""
import requests
from Bio import SeqIO
from io import StringIO
import sys

gene = sys.argv[1]
taxon = sys.argv[2]
output = sys.argv[3]

query = f"gene:{gene} AND taxonomy_id:{taxon}"

url = "https://rest.uniprot.org/uniprotkb/search"

params = {
    "query": query,
    "format": "fasta",
    "size": 500
}

response = requests.get(url, params=params)

if response.status_code != 200:
    raise Exception("UniProt request failed")

fasta_data = response.text

records = list(SeqIO.parse(StringIO(fasta_data), "fasta"))

SeqIO.write(records, output, "fasta")

print(f"Downloaded {len(records)} sequences")
