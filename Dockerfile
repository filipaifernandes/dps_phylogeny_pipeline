FROM snakemake/snakemake:latest

RUN conda install -y -c conda-forge -c bioconda \
    biopython \
    requests \
    mafft \
    cd-hit \
    trimal \
    iqtree \
    && conda clean -a -y

ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /pipeline
