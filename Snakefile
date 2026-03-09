configfile: "config.yaml"

# List of proteins from config.yaml
PROTEINS = config["proteins"]

#### Final target of the workflow ####
rule all:
    input:
        "data/trees/final.treefile"

#### Fetch protein sequences from NCBI ####
rule fetch_ncbi:
    output:
        "data/raw/{protein}/ncbi.fasta"
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/fetch_sequences_ncbi.py \
        {wildcards.protein} \
        {config[taxon]} \
        {config[email]} \
        {config[retmax]} \
        {output}
        """

#### Fetch protein sequences from UniProt ####
rule fetch_uniprot:
    output:
        "data/raw/{protein}/uniprot.fasta"
    container:
        "docker://filipafernandes/dps_pipeline:005"
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/fetch_sequences_uniprot.py \
        {wildcards.protein} \
        {config[taxon]} \
        {config[email]} \
        {config[retmax]} \
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
        python scripts/merge_clean_fasta.py {input[0]} {input[1]} {output}
        """

rule combine_proteins:
    input:
        expand("data/cleaned/{protein}/cleaned.fasta", protein=PROTEINS)
    output:
        "data/combined/all_sequences.fasta"
    shell:
        """
        mkdir -p data/combined
        cat {input} > {output}
        """

#### Reduce redundancy - using CD-HIT ####
rule cdhit:
    input:
        "data/cleaned/{protein}/cleaned.fasta"
    output:
        "data/cleaned/{protein}/nonredundant.fasta"
    shell:
        """
        mkdir -p $(dirname {output})
        cd-hit -i {input} -o {output} -c {config["cdhit_identity"]} -n 5
        """


#### Multiple sequence alignment - using MAFFT ####
rule align:
    input:
        "data/combined/all_sequences.fasta"
    output:
        "data/aligned/aligned.fasta"
    conda:
        "envs/mafft.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        # limpa os cabeçalhos antes de alinhar
        awk '/^>/ {{print $1}} /^[^>]/ {{print $0}}' {input} > data/combined/all_sequences_clean.fasta
        mafft --linsi data/combined/all_sequences_clean.fasta > {output}
        """

#### Alignment trimming - using TrimAl ####
rule trim:
    input:
        "data/aligned/aligned.fasta"
    output:
        "data/aligned/aligned_trimmed.fasta"
    conda:
        "envs/trimal.yaml"
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """

#### Phylogenetic inference - using IQ-TREE ####
rule tree:
    input:
        "data/aligned/aligned.fasta"
    output:
        "data/trees/final.treefile"
    conda:
        "envs/iqtree.yaml"
    shell:
        """
        mkdir -p data/trees
        iqtree2 -s {input} -m MFP \
	-bb {config[iqtree][bootstrap]} \
	-alrt {config[iqtree][alrt]} \
	-nt {config[iqtree][threads]} -redo
        mv data/aligned/aligned.fasta.treefile {output}
        """
