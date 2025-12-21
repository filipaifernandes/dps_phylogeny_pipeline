configfile: "config.yaml"

# List of proteins from config.yaml
PROTEINS = config["proteins"]

#### Final target of the workflow ####
rule all:
    input:
        expand("data/trees/{protein}.treefile", protein=PROTEINS)


#### Fetch protein sequences from NCBI ####
rule fetch_ncbi:
    output:
        "data/raw/{protein}/ncbi.fasta"
    conda:
        "envs/biopython.yaml"
    run:
        protein = wildcards.protein
        taxon = config["taxon"]
        email = config["email"]
        max_seqs = config["retmax"]
        out_file = output[0]

        shell(f"mkdir -p $(dirname {out_file}) && "
              f"python scripts/fetch_sequences_ncbi.py {protein} {taxon} {email} {max_seqs} {out_file}")


#### Fetch protein sequences from UniProt ####
rule fetch_uniprot:
    output:
        "data/raw/{protein}/uniprot.fasta"
    conda:
        "envs/requests.yaml"
    run:
        protein = wildcards.protein
        taxon = config["taxon"]
        out_file = output[0]

        shell(f"mkdir -p $(dirname {out_file}) && "
              f"python scripts/fetch_sequences_uniprot.py {protein} {taxon} {out_file}")


#### Merge and clean sequences ####
rule merge_clean:
    input:
        "data/raw/{protein}/ncbi.fasta",
        "data/raw/{protein}/uniprot.fasta"
    output:
        "data/cleaned/{protein}/cleaned.fasta"
    conda:
        "envs/biopython.yaml"
    run:
        inp = " ".join(input)
        out_file = output[0]

        shell(f"mkdir -p $(dirname {out_file}) && "
              f"python scripts/merge_and_clean_fasta.py {inp} {out_file}")


#### Reduce redundancy - using CD-HIT ####
rule cdhit:
    input:
        "data/cleaned/{protein}/cleaned.fasta"
    output:
        "data/cleaned/{protein}/nonredundant.fasta"
    conda:
        "envs/cdhit.yaml"
    run:
        inp = input[0]
        out_file = output[0]
        identity = config["cdhit_identity"]

        shell(f"mkdir -p $(dirname {out_file}) && "
              f"cd-hit -i {inp} -o {out_file} -c {identity} -n 5")


#### Multiple sequence alignment - using MAFFT ####
rule align:
    input:
        "data/cleaned/{protein}/nonredundant.fasta"
    output:
        "data/aligned/{protein}/aligned.fasta"
    conda:
        "envs/mafft.yaml"
    run:
        inp = input[0]
        out_file = output[0]

        shell(f"mkdir -p $(dirname {out_file}) && mafft --auto {inp} > {out_file}")


#### Alignment trimming - using TrimAl ####
rule trim:
    input:
        "data/aligned/{protein}/aligned.fasta"
    output:
        "data/aligned/{protein}/trimmed.fasta"
    conda:
        "envs/trimal.yaml"
    run:
        inp = input[0]
        out_file = output[0]

        shell(f"mkdir -p $(dirname {out_file}) && trimal -automated1 -in {inp} -out {out_file}")


#### Phylogenetic inference - using IQ-TREE ####
rule iqtree:
    input:
        "data/aligned/{protein}/trimmed.fasta"
    output:
        "data/trees/{protein}.treefile"
    conda:
        "envs/iqtree.yaml"
    run:
        inp = input[0]
        out_file = output[0]
        bootstrap = config["iqtree"]["bootstrap"]
        alrt = config["iqtree"]["alrt"]

        shell(f"mkdir -p $(dirname {out_file}) && "
              f"iqtree2 -s {inp} -m MFP -bb {bootstrap} -alrt {alrt} -nt AUTO")