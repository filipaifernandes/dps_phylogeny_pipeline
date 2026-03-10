configfile: "config.yaml"

# List of proteins from config.yaml
PROTEINS = config["proteins"]

#### Final target of the workflow ####
rule all:
    input:
        "data/trees/final.treefile",
        expand("data/trees/{protein}.treefile", protein=PROTEINS)


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
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/merge_clean_fasta.py {input[0]} {input[1]} {output}
        """


#### Reduce redundancy - using CD-HIT ####
rule cdhit:
    input:
        "data/cleaned/{protein}/cleaned.fasta"
    output:
        "data/cleaned/{protein}/nonredundant.fasta"
    conda:
        "envs/pipeline.yaml"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output})
        cd-hit -i {input} -o {output} -c {config[cdhit_identity]} -n {config[word_length]}
        """


#### Combine proteins ####
rule combine_proteins:
    input:
        expand("data/cleaned/{protein}/nonredundant.fasta", protein=PROTEINS)
    output:
        "data/combined/all_sequences.fasta"
    shell:
        """
        mkdir -p data/combined
        cat {input} > {output}
        """

#### Multiple sequence alignment - using MAFFT ####
rule align_combined:
    input:
        "data/combined/all_sequences.fasta"
    output:
        "data/aligned/aligned.fasta"
    conda:
        "envs/pipeline.yaml"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output})
        mafft --auto --thread {threads} {input} > {output}
        """


rule align_individual:
    input:
        "data/cleaned/{protein}/nonredundant.fasta"
    output:
        "data/aligned/{protein}_aligned.fasta"
    conda:
        "envs/pipeline.yaml"
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output})
        mafft --auto --thread {threads} {input} > {output}
        """


#### Alignment trimming - using TrimAl ####
rule trim_combined:
    input:
        "data/aligned/aligned.fasta"
    output:
        "data/aligned/aligned_trimmed.fasta"
    conda:
        "envs/pipeline.yaml"
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """


rule trim_individual:
    input:
        "data/aligned/{protein}_aligned.fasta"
    output:
        "data/aligned/{protein}_aligned_trimmed.fasta"
    conda:
        "envs/pipeline.yaml"
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """


#### Phylogenetic inference - using IQ-TREE ####
rule iqtree_combined:
    input:
        "data/aligned/aligned_trimmed.fasta"
    output:
        "data/trees/final.treefile"
    threads: 4
    conda:
        "envs/pipeline.yaml"
    shell:
        """
        mkdir -p data/trees
        iqtree2 -s {input} -m MFP \
        -bb {config[bootstrap]} \
        -alrt {config[alrt]} \
        -seed {config[seed]} \
        -nt {threads} -redo
        mv data/aligned/aligned_trimmed.fasta.treefile {output}
        """


rule iqtree_individual:
    input:
        "data/aligned/{protein}_aligned_trimmed.fasta"
    output:
        "data/trees/{protein}.treefile"
    threads: 4
    conda:
        "envs/pipeline.yaml"
    shell:
        """
        mkdir -p data/trees
        iqtree2 -s {input} -m MFP \
        -bb {config[iqtree][bootstrap]} \
        -alrt {config[iqtree][alrt]} \
        -nt {threads} -redo
        mv {input}.treefile {output}
        """
