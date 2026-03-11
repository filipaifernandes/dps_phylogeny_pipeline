FROM snakemake/snakemake:latest

RUN conda install -y -c conda-forge -c bioconda \
    biopython \
    requests \
    mafft \
    cd-hit \
    trimal \
    iqtree \
    seqkit \
    && conda clean -a -y

WORKDIR /pipeline
