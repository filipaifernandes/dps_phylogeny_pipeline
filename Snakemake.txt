configfile: "config.yaml"
PROTEINS = config["proteins"]

rule all:
    input:
        expand("data/trees/{protein}.treefile", protein=PROTEINS),
        expand("results/{protein}_alignment_stats.txt", protein=PROTEINS),
        expand("results/{protein}.lmap", protein=PROTEINS),
        expand("results/{protein}_blast_hits.tsv", protein=PROTEINS),
        expand("data/domains/{protein}.pfam.tbl", protein=PROTEINS)

#Fetch NCBI sequences
rule fetch_ncbi:
    output: "data/raw/{protein}_ncbi.fasta"
    conda: "envs/biopython.yaml"
    params: taxon=config["taxon"], email=config["email"], retmax=config["retmax"]
    shell: "python scripts/fetch_sequences_ncbi.py {wildcards.protein} {params.taxon} {params.email} {params.retmax} {output}"

#Fetch UniProt sequences
rule fetch_uniprot:
    output: "data/raw/{protein}_uniprot.fasta"
    conda: "envs/requests.yaml"
    params: taxon=config["taxon"]
    shell: "python scripts/fetch_sequences_uniprot.py {wildcards.protein} {params.taxon} {output}"

#Merge sequences from all databases
rule merge_sequences:
    input: 
        ncbi="data/raw/{protein}_ncbi.fasta",
        uniprot="data/raw/{protein}_uniprot.fasta"
    output: "data/raw/{protein}_combined.fasta"
    conda: "envs/biopython.yaml"
    shell: "python scripts/merge_sequences.py {input.ncbi} {input.uniprot} {output}"

#Clean merged sequences
rule clean_sequences:
    input: "data/raw/{protein}_combined.fasta"
    output: "data/cleaned/{protein}.clean.fasta"
    conda: "envs/biopython.yaml"
    params: min_len=config["min_length"], max_len=config["max_length"]
    shell: "python scripts/clean_sequences.py {input} {output} {params.min_len} {params.max_len}"

#Redundancy reduction
rule cdhit:
    input: "data/cleaned/{protein}.clean.fasta"
    output: "data/cleaned/{protein}.nr.fasta"
    conda: "envs/cdhit.yaml"
    params: c=config["cdhit_identity"]
    shell: "cd-hit -i {input} -o {output} -c {params.c}"

#Domain validation
rule hmmer_scan:
    input: "data/cleaned/{protein}.nr.fasta"
    output: "data/domains/{protein}.pfam.tbl"
    conda: "envs/hmmer.yaml"
    params: pfam="Pfam-A.hmm"
    shell: "hmmscan --tblout {output} {params.pfam} {input}"

#Alignment
rule align:
    input: "data/cleaned/{protein}.nr.fasta"
    output: "data/aligned/{protein}.aln.fasta"
    conda: "envs/mafft.yaml"
    shell: "mafft --maxiterate 1000 --localpair {input} > {output}"

#Alignment stats
rule alignment_stats:
    input: "data/aligned/{protein}.aln.fasta"
    output: "results/{protein}_alignment_stats.txt"
    conda: "envs/amas.yaml"
    shell: "amas summary -f fasta -i {input} > {output}"

#Trimming
rule trim:
    input: "data/aligned/{protein}.aln.fasta"
    output: "data/trimmed/{protein}.trim.fasta"
    conda: "envs/trimal.yaml"
    shell: "trimal -in {input} -out {output} -automated1"

#Phylogenetic tree
rule iqtree:
    input: "data/trimmed/{protein}.trim.fasta"
    output: "data/trees/{protein}.treefile"
    conda: "envs/iqtree.yaml"
    params: bb=config["iqtree"]["bootstrap"], alrt=config["iqtree"]["alrt"], nt=config["iqtree"]["threads"]
    shell: "iqtree2 -s {input} -m MFP -bb {params.bb} -alrt {params.alrt} -nt {params.nt}"

#Likelihood mapping
rule likelihood_mapping:
    input: "data/trimmed/{protein}.trim.fasta"
    output: "results/{protein}.lmap"
    conda: "envs/iqtree.yaml"
    shell: "iqtree2 -s {input} -lmap 10000"

#BLAST orthology check
rule blast_orthology_check:
    input: "data/cleaned/{protein}.nr.fasta"
    output: "results/{protein}_blast_hits.tsv"
    conda: "envs/blast.yaml"
    params: db="blast_db/Dps_refs"
    shell: "blastp -query {input} -db {params.db} -outfmt 6 -max_target_seqs 1 > {output}"
