FROM snakemake/snakemake:latest

RUN mamba install -y -c conda-forge -c bioconda \
    biopython \
    requests \
    mafft \
    cd-hit \
    trimal \
    iqtree \
    && mamba clean -a -y

WORKDIR /pipeline
