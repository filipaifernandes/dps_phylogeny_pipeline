configfile: "config.yaml"
container: "docker://filipafernandes/dps_pipeline:005"

PROTEINS = config["proteins"]
TRUNCATIONS = config["truncations"]


rule all:
    input:
        expand("data/raw/{protein}/ncbi.fasta", protein=PROTEINS),
        expand("data/raw/{protein}/uniprot.fasta", protein=PROTEINS),
        expand("data/cleaned/{protein}/{protein}_trunc.fasta", protein=TRUNCATIONS.keys()),
        expand("data/cleaned/{protein}/nonredundant.fasta", protein=PROTEINS),
        expand("data/aligned/{protein}_aligned.fasta", protein=PROTEINS),
        expand("data/aligned/{protein}_aligned_trimmed.fasta", protein=PROTEINS),
        expand("data/trees/{protein}.treefile", protein=PROTEINS),
        "data/trees/final_trunc.treefile",
        "data/trees/final_full.treefile"

#Fetch NCBI
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

#Fetch UniProt
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

#Merge and clean
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

#CD-HIT
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

#Truncation
rule truncate:
    input:
        "data/cleaned/{protein}/nonredundant.fasta"
    output:
        "data/cleaned/{protein}/{protein}_trunc.fasta"
    params:
        start = lambda wc: config["truncations"][wc.protein]["start"],
        end   = lambda wc: config["truncations"][wc.protein]["end"]
    shell:
        """
        python scripts/truncate_protein.py {input} {output} {params.start} {params.end}
        """

#Combine: trunc + full
rule combine_trunc_full:
    input:
        expand("data/cleaned/{protein}/nonredundant.fasta", protein=PROTEINS),
        "data/cleaned/dps1/dps1_trunc.fasta"
    output:
        "data/combined/all_sequences_trunc.fasta"
    shell:
        """
        mkdir -p data/combined
        cat {input} > {output}
        """

#Combine: full only
rule combine_full_only:
    input:
        expand("data/cleaned/{protein}/nonredundant.fasta", protein=PROTEINS)
    output:
        "data/combined/all_sequences_full.fasta"
    shell:
        """
        mkdir -p data/combined
        cat {input} > {output}
        """

#Alignment
rule align_trunc_full:
    input:
        "data/combined/all_sequences_trunc.fasta"
    output:
        "data/aligned/aligned_trunc.fasta"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output})
        mafft {config[mafft][method]} {input} > {output}
        """

rule align_full:
    input:
        "data/combined/all_sequences_full.fasta"
    output:
        "data/aligned/aligned_full.fasta"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output})
        mafft {config[mafft][method]} {input} > {output}
        """

#Trimming
rule trim_trunc:
    input:
        "data/aligned/aligned_trunc.fasta"
    output:
        "data/aligned/aligned_trunc_trimmed.fasta"
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """

rule trim_full:
    input:
        "data/aligned/aligned_full.fasta"
    output:
        "data/aligned/aligned_full_trimmed.fasta"
    shell:
        """
        trimal -in {input} -out {output} -automated1
        """

#Trees
rule iqtree_trunc:
    input:
        "data/aligned/aligned_trunc_trimmed.fasta"
    output:
        "data/trees/final_trunc.treefile"
    threads: 4
    shell:
        """
        iqtree -s {input} \
        -m MFP \
        -B {config[iqtree][bootstrap]} \
        --alrt {config[iqtree][alrt]} \
        -T AUTO \
        --prefix data/trees/final_trunc \
        --redo
        """

rule iqtree_full:
    input:
        "data/aligned/aligned_full_trimmed.fasta"
    output:
        "data/trees/final_full.treefile"
    threads: 4
    shell:
        """
        iqtree -s {input} \
        -m MFP \
        -B {config[iqtree][bootstrap]} \
        --alrt {config[iqtree][alrt]} \
        -T AUTO \
        --prefix data/trees/final_full \
        --redo
        """

#Individual trees

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

rule trim_individual:
    input:
        "data/aligned/{protein}_aligned.fasta"
    output:
        "data/aligned/{protein}_aligned_trimmed.fasta"
    shell:
        """
        trimal -in {input} -out {output} -automated1
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
        -T AUTO \
        --prefix data/trees/{wildcards.protein} \
        --redo
        """
