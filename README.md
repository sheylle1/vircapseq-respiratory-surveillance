# Implementation of VirCapSeq-VERT Metagenomic Surveillance for Respiratory Viruses in South Africa

This repository contains the code and de-identified metadata used for the analyses presented in the manuscript:

> **Implementation of hybridisation capture metagenomic surveillance for respiratory viruses in South Africa**

## Overview

This study evaluates the implementation of the VirCapSeq-VERT hybridisation capture metagenomic sequencing assay for respiratory virus surveillance in South Africa.The repository supports reproducible research.

*Note*: Materials in this repository may be updated upon final publication.


## Repository Contents

```
scripts/          Analysis scripts
data/             De-identified sample metadata
results/          Generated figures and summary tables
```

## Data Availability

Human-read-removed FASTQ files are publicly available from the NCBI Sequence Read Archive (SRA):

**BioProject:** PRJNA1455986

Consensus genome sequences have been submitted to GISAID (accession numbers provided in Supplementary Table X).
Only de-identified metadata required to reproduce the analyses is included in this repository.

## Reproducing the Analysis

1. Download the human-read-removed FASTQ files from the SRA.
2. Obtain the consensus genome sequences from GISAID (where permitted).
3. Clone this repository.
4. Install the required software and dependencies.
5. Run the analysis scripts in the order indicated in the `scripts/` directory.

## Software

Analyses were performed using:

* R (version X.X.X)
* Python (version X.X)
* Genome Detective (v2.17–2.21)
* Nextclade (v3.21.0)
* FastQC (v0.11.2)
* MultiQC

Additional package versions are provided in the individual analysis scripts.


## Contact

For questions or issues regarding this repository, please contact the corresponding author listed in the manuscript.
