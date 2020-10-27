# Pan-cancer analysis of genomic profiles and immune landscape of Pembrolizumab-treated metastatic solid tumors.

## Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Demo](#demo)
- [Results](#results)
- [License](./LICENSE)
- [Issues](https://github.com/pughlab/inspire-genomics/issues)
- [Citation](#citation)

# Overview
Pan-cancer assessments of predictive features and modulations in the tumor and immune microenvironment of patients undergoing immune checkpoint blockade are limited. Atypical clinical responses restrict accurate patient stratification for biomarker discovery and evaluation. Here, we present an in-depth analysis of clinical, whole-exome, transcriptome, and circulating tumor DNA (ctDNA) profiles of 73 patients with advanced solid tumors, across 30 cancer types, from a phase II clinical trial of pembrolizumab (NCT02644369). Patients stratified by ctDNA and tumor burden dynamics corresponded with survival and clinical benefit (CB). High mutation burden, high expression of immune signatures, and mutations in BRCA2 were associated with pembrolizumab molecular sensitivity, while abundant copy-number alterations and B2M loss-of-heterozygosity corresponded with resistance. Upon treatment, induction of genes expressed by T-, B-, and myeloid cell populations were consistent with sensitivity and resistance. We identified PLA2G2D, an immune-regulating phospholipase, as a potential therapeutic target and marker of immune resistance. Together, these findings provide insights into the diversity of immunogenomic mechanisms that underpin pembrolizumab outcomes.

# Repo Contents
This repo contains custom R scripts to reproduce the main analyses and figures presented in the submitted manuscript to Nature Communications. Analyses include: 
- Waterfall plots
- Survival analysis
- Co-mut plot
- Mutation enrichment
- ssGSEA (single sample geneset enrichment)
- Differential gene-expression

# System Requirements
## Hardware Requirements
All scripts in this repository was developed on a system with the following specs:

RAM: 16 GB
CPU: 2.6 GHz Intel Core i7
MAC OSX: 10.9.5

## Software Requirements
General software: R (version 3.3.1)

### R-packages
Visualization: pheatmap, maftools 
RNAseq analysis: RSEM (version 1.3.0), sva (version 3.36), DESeq2 (version 3.11), CIBERSORT (https://cibersort.stanford.edu/), gsva(veresion 3.11),

# Demo
# Results

# Citation
For usage of the package and associated manuscript, please cite according to the enclosed [citation.bib](./citation.bib).
