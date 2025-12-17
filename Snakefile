configfile: "config.yaml"

#Final target of the workflow
rule all:
    input:
        "data/trees/dps.treefile"


#### Fetch protein sequences from NCBI ####
rule fetch_ncbi:
    output:
        "data/raw/ncbi.fasta"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/fetch_sequences_ncbi.py \
            {config[gene]} \
            {config[taxon]} \
            {config[email]} \
            {config[max_seqs]} \
            {output}
        """


#### Fetch protein sequences from UniProt ####
rule fetch_uniprot:
    output:
        "data/raw/uniprot.fasta"
    conda:
        "envs/requests.yaml"
    shell:
        """
        python scripts/fetch_sequences_uniprot.py \
            {config[gene]} \
            {config[taxon]} \
            {output}
        """


#### Merge and clean sequences ####
rule merge_clean:
    input:
        "data/raw/ncbi.fasta",
        "data/raw/uniprot.fasta"
    output:
        "data/cleaned/cleaned.fasta"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/merge_and_clean_fasta.py \
            {input} \
            {output}
        """


#### Reduce redundancy - using CD-HIT
rule cdhit:
    input:
        "data/cleaned/cleaned.fasta"
    output:
        "data/cleaned/nonredundant.fasta"
    conda:
        "envs/cdhit.yaml"
    shell:
        """
        cd-hit -i {input} -o {output} -c 0.95 -n 5
        """

#### Multiple sequence alignment ####
rule align:
    input:
        "data/cleaned/nonredundant.fasta"
    output:
        "data/aligned/aligned.fasta"
    conda:
        "envs/mafft.yaml"
    shell:
        """
        mafft --auto {input} > {output}
        """


#### Alignment trimming - using TrimAl ####
rule trim:
    input:
        "data/aligned/aligned.fasta"
    output:
        "data/aligned/trimmed.fasta"
    conda:
        "envs/trimal.yaml"
    shell:
        """
        trimal -automated1 -in {input} -out {output}
        """


#### Phylogenetic inference - using IQ-TREE ####
rule iqtree:
    input:
        "data/aligned/trimmed.fasta"
    output:
        "data/trees/dps.treefile"
    conda:
        "envs/iqtree.yaml"
    shell:
        """
        iqtree2 -s {input} -m MFP -bb 1000 -alrt 1000
        """