# DPS Phylogeny Pipeline

**Automated and reproducible workflow for the retrieval, curation, alignment, and phylogenetic analysis of Dps proteins in Deinococcus species.**

This pipeline performs a complete bioinformatics workflow to analyze Dps1 and Dps2 proteins from Deinococcus species, with a particular interest in Deinococcus radiodurans.
The workflow is made by using Snakemake and executed within a Docker container, ensuring a fully reproducible and portable computational environment.

The pipeline automatically retrieves sequences from multiple databases, performs dataset preprocessing, generates a multiple sequence alignment, and reconstructs a maximum-likelihood phylogenetic tree.

---

## Table of Contents

1. [Overview](#overview)  
2. [Key Features](#key-features)  
3. [Installation](#installation)  
4. [Usage](#usage)  
5. [Pipeline Workflow](#pipeline-workflow)  
6. [Reproducibility and Automation](#reproducibility-and-automation)  
7. [Continuous Integration](#continuous-integration)  
8. [Docker Execution](#optional-docker-execution)  
9. [Configuration](#configuration)  
10. [Directory Structure](#directory-structure)  
11. [References](#references)  
12. [Contact](#contact)  

---

## Overview

The DPS Phylogeny Pipeline is designed to fully provide a **reproducible and automated workflow for phylogenetic analysis of Dps proteins**. 
The pipeline performs the following steps:

- Sequence retrieval from **NCBI** and **UniProt**  
- Sequence merging and cleaning, including duplicate's removal  
- Combination of Dps1 and Dps2 datasets
- Redundancy reduction using **CD-HIT**  
- Multiple sequence alignment using **MAFFT**  
- Alignment trimming using **TrimAl**  
- Maximum-likelihood phylogenetic inference with **IQ-TREE**, including model selection and branch support

All steps are executed automatically through Snakemake rules, ensuring that **no manual intervention is required**.

---

## Key Features

- **Automated workflow:** Single-command execution of the complete pipeline  
- **Multi-database retrieval:** Fetches sequences from NCBI and UniProt  
- **High-quality preprocessing:** duplication, filtering, trimming  
- **Protein dataset integration:** combines multiple Dps proteins into a unified dataset
- **High-quality alignment:** generated with MAFFT
- **Phylogenetic rigor:** Alignment trimming, model selection, bootstrap and aLRT support  
- **Containerized execution:** Docker ensures identical behaviour across systems  
- **CI validation:** Lightweight GitHub Actions workflow ensures workflow integrity on every commit  

---

## Installation
The pipeline uses **Snakemake with Singularity containers** to guarantee reproducibility.

1. Install **Snakemake**
Recommended installation via Conda:
```bash
conda install -c conda-forge -c bioconda snakemake
```

2. Install **Singularity**
```bash
sudo apt update
sudo apt install singularity-container
```

3. Clone the repository:

```bash
git clone https://github.com/yourname/dps_phylogeny_pipeline.git
cd dps_phylogeny_pipeline
```

4. Run the pipeline:

```bash
snakemake --use-singularity --cores 4
```
The required container image 
```(docker://filipafernandes/dps_pipeline:005)
```
 will be automatically downloaded and executed through Singularity.
---

## Usage

- Only one command is required for complete execution.
- Outputs are structured for clarity:

```
data/
├── raw/             # Retrieved sequences
├── cleaned/         # Merged and filtered sequences
├── aligned/         # MAFFT alignments and trimmed alignments
└── trees/           # Phylogenetic trees
```

- The pipeline can be run multiple times.

---

## Pipeline Workflow

### 1. Sequence Retrieval
Protein sequences are retrieved from two sources:

**NCBI** 
- Retrieves usinf the **Biopython Entrez API**  
- Queries constructed using the configured protein names and taxon

**UniProt:** 
- Retrieves using the **UniProt REST API**  
- Complementary dataset to increase sequence coverage

Sequences are stored in:
data/raw/{protein}/

### 2. Sequence Cleaning
Sequences from both databases are merged and processed using a Python script.
The cleaning step currently performs:
- Merge fasta files form NCBI and UniProt  
- Remove duplicates and ambiguous residues  
- Generate a clean FASTA file for downstream analysis

The resulting dataset is stored in:
data/cleaned/{protein}/cleaned.fasta

Cleaned sequences for each protein are combined into a single dataset:
data/combined/all_sequences.fasta
This allows phylogenetic analysis across multiple Dps homologs.

### 3. Redundancy Reduction

- **CD-HIT** collapses sequences with ≥95% identity  
- Reduces bias from overrepresented sequences

### 4. Multiple Sequence Alignment
Sequences are aligned using MAFFT.

The workflow uses the MAFFT automatic mode:
```
mafft --auto
```

This allows MAFFT to automatically select the most appropriate alignment algorithm based on dataset characteristics.

The resulting alignment is stored in:
data/aligned/aligned.fasta

### 5. Alignment Trimming
The alignment is processed using TrimAl to remove poorly aligned regions.

The pipeline uses:
```
trimal -automated1
```

This mode automatically determines trimming thresholds to improve alignment quality.

Output:
data/aligned/aligned_trimmed.fasta

### 6. Phylogenetic Inference
Phylogenetic reconstruction is performed using IQ-TREE under a maximum-likelihood framework.

The pipeline performs:
- ModelFinder Plus (MFP) for substitution model selection
- Ultrafast bootstrap analysis
- SH-aLRT branch support estimation

Example command used by the workflow:
```
iqtree2 -s alignment.fasta -m MFP -bb 1000 -alrt 1000
```

The resulting tree is saved as:
data/trees/final.treefile

---

## Reproducibility and Automation
Reproducibility is ensured through:

- Snakemake workflow management
- Docker containerized environment
- Clearly defined configuration parameters
- Automated dependency resolution

Snakemake also constructs a Directed Acyclic Graph (DAG) of the workflow to guarantee correct execution order.

---

## Continuous Integration
The repository includes a GitHub Actions workflow that automatically validates the pipeline.

The CI workflow performs:
- Snakemake dry-run execution
- Validation of the Snakefile
- Verification of rule dependencies and DAG structure

This ensures that changes to the repository do not break the workflow.

---

## Docker Execution
The pipeline uses a Docker container based on:
```
snakemake/snakemake:latest
```

Additional bioinformatics tools are installed in the image:
- Biopython
- Requests
- MAFFT
- CD-HIT
- TrimAl
- IQ-TREE

This ensures the pipeline runs identically on Linux, macOS, or Windows systems.

---

## Configuration
Pipeline behaviour is controlled through the `config.yaml` file.

Example configuration:
```yaml
proteins: 
- Dps1 
- Dps2 

taxon: Deinococcus 
focus_species: Deinococcus radiodurans 

min_length: 150 
max_length: 300 

cdhit_identity: 0.95 

email: youremail@example.pt 
retmax: 500
```

**Parameter description**
| Parameter      | Description                                       |
| -------------- | ------------------------------------------------- |
| proteins       | Target proteins to retrieve                       |
| taxon          | Taxonomic group used for sequence search          |
| focus_species  | Species of special interest                       |
| min_length     | Minimum allowed sequence length                   |
| max_length     | Maximum allowed sequence length                   |
| cdhit_identity | Sequence identity threshold for CD-HIT clustering |
| email          | Required for NCBI Entrez API                      |
| retmax         | Maximum number of sequences retrieved per query   |

---

## Directory Structure

The workflow structure generated by Snakemake is shown below:

![Pipeline DAG](/home/pipas/dps_phylogeny_pipeline/dag.png)

---

## References

- **Snakemake**: Köster & Rahmann, Bioinformatics 2012  
- **Biopython**: Cock et al., Bioinformatics 2009  
- **MAFFT**: Katoh & Standley, Mol Biol Evol 2013  
- **TrimAl**: Capella-Gutiérrez et al., Bioinformatics 2009  
- **IQ-TREE**: Minh et al., Mol Biol Evol 2020  
- **CD-HIT**: Li & Godzik, Bioinformatics 2006  

---

## Contact

Filipa Fernandes

Bioinformatics Student

For questions, troubleshooting or contributions regarding the pipeline:
[filipaifernandes.2005@gmail.com](mailto:filipaifernandes.2005@gmail.com)

