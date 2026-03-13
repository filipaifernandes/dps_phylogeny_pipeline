configfile: "config.yaml"
container: "docker://filipafernandes/dps_pipeline:005"
PROTEINS = config["proteins"]



#### Final target of the workflow ####
rule all:
    input:
        expand("data/raw/{protein}/ncbi.fasta", protein=PROTEINS),
        expand("data/raw/{protein}/uniprot.fasta", protein=PROTEINS),
        "data/cleaned/dps1/dps1_trunc.fasta",
        expand("data/cleaned/{protein}/nonredundant.fasta", protein=PROTEINS),
        "data/aligned/aligned.fasta",
        expand("data/aligned/{protein}_aligned.fasta", protein=PROTEINS),
        "data/aligned/aligned_trimmed.fasta",
        expand("data/aligned/{protein}_aligned_trimmed.fasta", protein=PROTEINS),
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
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output})
        cd-hit -i {input} -o {output} -c {config[cdhit_identity]} -n {config[word_length]}
        """

##Exclusive to this analysis (can be documented if not needed) 
#### Create truncated dps1 sequence (aa 54–207) ####
rule truncate_dps1:
    input:
        "data/cleaned/dps1/nonredundant.fasta"
    output:
        "data/cleaned/dps1/dps1_trunc.fasta"
    shell:
        """
        python scripts/truncate_dps1.py {input} {output} {config[truncation][start]} {config[truncation][end]}
        """

#### Combine proteins ####
rule combine_proteins:
    input:
        expand("data/cleaned/{protein}/nonredundant.fasta", protein=PROTEINS),
        "data/cleaned/dps1/dps1_trunc.fasta"
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
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output})
	mafft --{config[mafft][method]}
        """


rule align_individual:
    input:
        "data/cleaned/{protein}/nonredundant.fasta"
    output:
        "data/aligned/{protein}_aligned.fasta"
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
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """


rule trim_individual:
    input:
        "data/aligned/{protein}_aligned.fasta"
    output:
        "data/aligned/{protein}_aligned_trimmed.fasta"
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
    shell:
        """
        iqtree -s {input} \
        -m MFP \
        -B {config[iqtree][bootstrap]} \
        --alrt {config[iqtree][alrt]} \
        -T {threads} \
        --prefix data/trees/final \
        --redo
        """

rule iqtree_individual:
    input:
        "data/aligned/{protein}_aligned_trimmed.fasta"
    output:
        "data/trees/{protein}.treefile"
    threads: 4
    shell:
        """
        iqtree -s {input} \
        -m MFP \
        -B {config[iqtree][bootstrap]} \
        --alrt {config[iqtree][alrt]} \
        -T {threads} \
        --prefix data/trees/{wildcards.protein} \
        --redo
        """
