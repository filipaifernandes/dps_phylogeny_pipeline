# DPS Phylogeny Pipeline

**Automated, reproducible workflow for the retrieval, curation, alignment, and phylogenetic analysis of Dps proteins in Deinococcus species.**

This pipeline is designed to streamline bioinformatics analysis of **Dps1 and Dps2 proteins of Deinococcus species**, with a primary focus on *Deinococcus radiodurans*. Leveraging **Snakemake** for workflow management, **Conda** for environment reproducibility, and optional **Docker** for OS-level reproducibility, the pipeline integrates sequences from multiple sources, ensures high-quality data preprocessing, and performs robust phylogenetic inference.

---

## Table of Contents

1. [Overview](#overview)  
2. [Key Features](#key-features)  
3. [Installation](#installation)  
4. [Usage](#usage)  
5. [Pipeline Workflow](#pipeline-workflow)  
6. [Reproducibility and Automation](#reproducibility-and-automation)  
7. [Continuous Integration](#continuous-integration)  
8. [Optional Docker Execution](#optional-docker-execution)  
9. [Configuration](#configuration)  
10. [Directory Structure](#directory-structure)  
11. [References](#references)  
12. [Contact](#contact)  

---

## Overview

The DPS Phylogeny Pipeline is intended for **systematic and reproducible phylogenetic analysis** of Dps proteins. It automates the following steps:

- Sequence retrieval from **NCBI** and **UniProt**  
- Sequence merging and cleaning, including duplicate's removal  
- Redundancy reduction using **CD-HIT**  
- Multiple sequence alignment using **MAFFT**  
- Alignment trimming using **TrimAl**  
- Maximum-likelihood phylogenetic inference with **IQ-TREE**, including model selection and branch support

All steps are executed automatically using Snakemake rules, ensuring that **no manual intervention is required**.

---

## Key Features

- **Automated workflow:** Single-command execution of the complete pipeline  
- **Multi-database integration:** Fetches sequences from NCBI and UniProt  
- **High-quality preprocessing:** Deduplication, filtering, trimming  
- **Reproducible environments:** Each tool is installed in a dedicated Conda environment  
- **Phylogenetic rigor:** Alignment trimming, model selection, bootstrap and aLRT support  
- **Optional OS-level reproducibility:** Docker image for cross-platform consistency  
- **CI validation:** Lightweight GitHub Actions workflow ensures workflow integrity on every commit  

---

## Installation

### Native Installation

1. Install **Miniconda/Anaconda**.  
2. Clone the repository:

```bash
git clone https://github.com/yourname/dps_phylogeny_pipeline.git
cd dps_phylogeny_pipeline
```

3. Run the pipeline (tools installed automatically via Conda):

```bash
snakemake --use-conda --cores 4
```

### Optional Docker Execution

Docker ensures the pipeline runs identically on any system without manual setup:

```bash
docker build -t dps_pipeline .
docker run -it dps_pipeline
```

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

- The pipeline can be runed multiple times with **consistent results**.

---

## Pipeline Workflow

### 1. Sequence Retrieval

- **NCBI:** Biopython Entrez API  
- **UniProt:** REST API  
- Retrieves protein sequences matching gene and taxon parameters.

### 2. Sequence Cleaning

- Merge sequences from multiple sources  
- Remove duplicates and ambiguous residues  
- Generate a clean FASTA file for downstream analysis

### 3. Redundancy Reduction

- **CD-HIT** collapses sequences with ≥95% identity  
- Reduces bias from overrepresented sequences

### 4. Multiple Sequence Alignment

- **MAFFT** aligns non-redundant sequences  
- Ensures accurate residue positioning for phylogenetic inference

### 5. Alignment Trimming

- **TrimAl** removes poorly aligned regions  
- Improves signal-to-noise ratio for tree inference

### 6. Phylogenetic Inference

- **IQ-TREE** performs maximum-likelihood inference  
- **ModelFinder Plus (MFP)** selects the best-fit substitution model  
- **Bootstrap (1000 replicates) and aLRT** branch support values generated  

---

## Reproducibility and Automation

- Each bioinformatics tool resides in its **own Conda environment**  
- Optional Docker ensures **OS-level reproducibility**  
- Snakemake ensures **all steps are executed automatically** in correct order  
- Multiple runs produce **identical outputs** for the same input parameters

---

## Continuous Integration

- GitHub Actions workflow performs a **dry-run** on every commit and pull request  
- Validates:
  - Snakefile syntax  
  - Rule dependencies and DAG integrity  
  - Config file correctness  
- Ensures pipeline maintainability and error detection without performing heavy computations

---

## Optional Docker Execution

- Based on the **official Snakemake Docker image (`snakemake/snakemake:v7.32.4`)**  
- Provides an isolated execution environment  
- Ensures identical behaviour across machines, regardless of OS or installed packages

---

## Configuration

`config.yaml` allows customization of:

- `gene` — Gene name (Dps1 or Dps2)  
- `taxon` — Target organism (default: Deinococcus)  
- `email` — Required for NCBI Entrez API  
- `max_seqs` — Maximum sequences to retrieve  

Example:

```yaml
gene: "dps"
taxon: "Deinococcus"
email: "your_email@domain.com"
max_seqs: 500
```

---

## Directory Structure

```
dps_phylogeny_pipeline/
├── Snakefile                 # Main workflow
├── config.yaml               # User configuration
├── README.md                 # Documentation
├── Dockerfile                # Optional Docker execution
├── envs/                     # Conda environments per tool
├── scripts/                  # Python scripts for sequence retrieval and cleaning
├── data/                     # Raw, cleaned, aligned sequences and trees
└── .github/workflows/        # CI/CD configuration
```

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

For questions, troubleshooting, or contributions:  
**filipaifernandes.2005@gmail.com**
