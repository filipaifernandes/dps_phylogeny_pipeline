FROM snakemake/snakemake:v7.32.4

WORKDIR /pipeline

#Copy workflow
COPY . /pipeline

#Ensure conda is configured properly
RUN conda config --set channel_priority strict

#Run the pipeline
CMD ["snakemake", "--use-conda", "--cores", "4"]