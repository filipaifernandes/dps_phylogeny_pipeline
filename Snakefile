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
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/fetch_sequences_ncbi.py \
            {wildcards.protein} \
            {config["taxon"]} \
            {config["email"]} \
            {config["retmax"]} \
            {output}
        """


#### Fetch protein sequences from UniProt ####
rule fetch_uniprot:
    output:
        "data/raw/{protein}/uniprot.fasta"
    conda:
        "envs/requests.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/fetch_sequences_uniprot.py \
            {wildcards.protein} \
            {config['taxon']} \
            {output}
        """


#### Merge and clean sequences ####
rule merge_clean:
    input:
        "data/raw/{protein}/ncbi.fasta",
        "data/raw/{protein}/uniprot.fasta"
    output:
        "data/cleaned/{protein}/cleaned.fasta"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/merge_and_clean_fasta.py {input} {output}
        """


#### Reduce redundancy - using CD-HIT ####
rule cdhit:
    input:
        "data/cleaned/{protein}/cleaned.fasta"
    output:
        "data/cleaned/{protein}/nonredundant.fasta"
    conda:
        "envs/cdhit.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        cd-hit -i {input} -o {output} -c {config[cdhit_identity]} -n 5
        """


#### Multiple sequence alignment - using MAFFT ####
rule align:
    input:
        "data/cleaned/{protein}/nonredundant.fasta"
    output:
        "data/aligned/{protein}/aligned.fasta"
    conda:
        "envs/mafft.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        mafft --auto {input} > {output}
        """


#### Alignment trimming - using TrimAl ####
rule trim:
    input:
        "data/aligned/{protein}/aligned.fasta"
    output:
        "data/aligned/{protein}/trimmed.fasta"
    conda:
        "envs/trimal.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        trimal -automated1 -in {input} -out {output}
        """


#### Phylogenetic inference - using IQ-TREE ####
rule iqtree:
    input:
        "data/aligned/{protein}/trimmed.fasta"
    output:
        "data/trees/{protein}.treefile"
    conda:
        "envs/iqtree.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        iqtree2 -s {input} -m MFP -bb {config[iqtree][bootstrap]} -alrt {config[iqtree][alrt]} -nt AUTO
        """