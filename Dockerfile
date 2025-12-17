FROM snakemake/snakemake:v7.32.4

#Set working director
WORKDIR /pipeline

#Copy pipeline files into container
COPY . /pipeline

#Default command:
#uses Conda environments defined in the workflow
CMD ["snakemake", "--use-conda", "--cores", "4"]