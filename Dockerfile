FROM snakemake/snakemake:v7.32.4

# Install system tools needed
RUN apt-get update && apt-get install -y \
    mafft \
    iqtree \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workflow

COPY . /workflow

ENTRYPOINT ["snakemake"]
CMD ["--cores", "1"]
