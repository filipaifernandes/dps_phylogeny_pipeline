# DPS Phylogeny Pipeline

**Automated and reproducible workflow for the retrieval, curation, alignment, and phylogenetic analysis of Dps proteins in Deinococcus species.**

This pipeline performs a complete bioinformatics workflow to analyze Dps1 and Dps2 proteins from *Deinococcus* species, with a particular interest in *Deinococcus radiodurans*.
The workflow is made by using Snakemake and executed within a Docker container, ensuring a fully reproducible and portable computational environment.

The pipeline automatically retrieves sequences from multiple databases, performs dataset preprocessing, generates a multiple sequence alignment, and reconstructs a maximum-likelihood phylogenetic tree.

Two complementary phylogenetic analyses are performed:
- A **combined dataset including full-length and truncated Dps1 sequences**
- A **full-length-only dataset**

This enables direct comparison between conserved-domain-focused and full-protein evolutionary signals.

---

## Table of Contents

1. [Overview](#overview)
2. [Key Features](#key-features)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Output structure](#outputstructure)
6. [Pipeline Workflow](#pipeline-workflow)
7. [Reproducibility and Automation](#reproducibility-and-automation)
8. [Continuous Integration](#continuous-integration)
9. [Docker Image Construction](#docker-image-construction)
10. [Configuration](#configuration)
11. [Directory Structure](#directory-structure)
12. [License](#license)
13. [References](#references)
14. [Contact](#contact)

---

## Overview

The DPS Phylogeny Pipeline is designed to fully provide a **reproducible and automated workflow for phylogenetic analysis of Dps proteins**. 
The pipeline performs the following steps:

- Sequence retrieval from **NCBI** and **UniProt**
- Sequence merging and cleaning, including duplicate's removal
- Redundancy reduction using **CD-HIT**
- Dataset construction (full-length and truncated variants)
- Multiple sequence alignment using **MAFFT**
- Alignment trimming using **TrimAl**
- Maximum-likelihood phylogenetic inference with **IQ-TREE**, including model selection and branch support

All steps are executed automatically through Snakemake rules, ensuring that **no manual intervention is required**.

---

## Key Features

- **Automated workflow:** Single-command execution of the complete pipeline
- **Multi-database retrieval:** Fetches sequences from NCBI and UniProt
- **High-quality preprocessing:** duplication, filtering, trimming
- **Dual dataset strategy:** full-length vs truncated+full comparison
- **High-quality alignment:** generated with MAFFT
- **Phylogenetic rigor:** Alignment trimming, model selection, bootstrap and aLRT support
- **Containerized execution:** Docker ensures identical behaviour across systems
- **CI validation:** Lightweight GitHub Actions workflow ensures workflow integrity on every commit
- **Dps1 truncation:** automatically truncates Dps1 sequences (aa 54–207) for downstream analyses

---

## Installation
The pipeline uses **Snakemake with Singularity/Apptainer containers**.

1. Install **Snakemake**

Recommended installation via Conda:
```bash
conda install -c conda-forge -c bioconda snakemake
```

2. Install Apptainer

Download from:
```bash
https://github.com/apptainer/apptainer/releases/tag/v1.4.5
```
For Ubuntu, it's recommended getting the **apptainer_1.4.5_amd64.deb** package. <br>
Install it:
```bash
sudo apt install ./apptainer_1.4.5_amd64.deb
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
The required container image **(docker://filipafernandes/dps_pipeline:005)** will be automatically downloaded and executed through Singularity.

---

## Usage

- Single command execution
- Fully reproducible
- Safe to re-run (Snakemake handles dependencies)

---
## Output structure
```
data/
├── raw/        # Retrieved sequences
├── cleaned/    # Filtered and nonredundant sequences
├── combined/   # Dataset combinations
├── aligned/    # Alignments and trimmed alignments
└── trees/      # Phylogenetic trees
```

---

## Pipeline Workflow

### 1. Sequence Retrieval
Sequences are retrieved from:

**NCBI** 
- Biopython Entrez API
- Query based on protein name and taxon

**UniProt:** 
- REST API
- Expands dataset coverage

Output:<br>
data/raw/{protein}/

### 2. Sequence Cleaning

- Merge fasta files form NCBI and UniProt
- Remove duplicates and ambiguous residues
- Generate clean datasets

Output:<br>
data/cleaned/{protein}/cleaned.fasta

### 3. Redundancy Reduction

- Performed using CD-HIT
- Collapses sequences with ≥95% identity
- Reduces redundancy bias

Output:<br>
data/cleaned/{protein}/nonredundant.fasta

### 4. Dps1 Truncation

- Extracts conserved region (aa 54–207)
- Applied only to Dps1

Output:<br>
data/cleaned/dps1/dps1_trunc.fasta

### 5. Dataset Construction

Two datasets are generated:

**Truncated + Full dataset**
- All full-length sequences
- Includes truncated Dps1

Output:<br>
data/combined/all_sequences_trunc.fasta

**Full-length dataset**

- Only full-length sequences
- data/combined/all_sequences_full.fasta

Output:<br>
data/combined/all_sequences_full.fasta

### 6. Multiple Sequence Alignment

Performed using MAFFT:
```
mafft --auto
```
Outputs:<br>
data/aligned/aligned_trunc.fasta
data/aligned/aligned_full.fasta

### 7. Alignment Trimming

Performed using TrimAl:

```
trimal -automated1
```

Outputs:
data/aligned/aligned_trunc_trimmed.fasta
data/aligned/aligned_full_trimmed.fasta

### 8. Phylogenetic Inference

Performed using IQ-TREE:

```
iqtree2 -s alignment.fasta -m MFP -bb 1000 -alrt 1000
```

Outputs:
**Truncated + Full tree**
data/trees/final_trunc.treefile

**Full-length tree**
data/trees/final_full.treefile

---
## Reproducibility and Automation
Reproducibility is ensured through:

- Snakemake workflow management
- Docker containerized execution
- Explicit configuration via `config.yaml`
- DAG-based execution ensures correct order

---

## Continuous Integration
GitHub Actions workflow performs:

- Snakemake dry-run
- DAG validation
- Rule consistency checks

---

## Docker Image Construction
Base image:

```
snakemake/snakemake:latest
```

Installed tools:
- Biopython
- Requests
- MAFFT
- CD-HIT
- TrimAl
- IQ-TREE

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

**Parameters**
| Parameter      | Description               |
| -------------- | ------------------------- |
| proteins       | Target proteins           |
| taxon          | Taxonomic group           |
| focus_species  | Species of interest       |
| min_length     | Minimum sequence length   |
| max_length     | Maximum sequence length   |
| cdhit_identity | CD-HIT identity threshold |
| email          | Required for NCBI API     |
| retmax         | Max sequences retrieved   |

---

## Directory Structure
The workflow structure generated by Snakemake is shown below:

![Pipeline DAG](dag.png)

Generate DAG:

```
snakemake --dag | dot -Tpng > dag.png
```

---

## License

This pipeline is licensed under the [MIT License](LICENSE).
You are free to use, modify and distribute it with proper attribution to the author.

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

