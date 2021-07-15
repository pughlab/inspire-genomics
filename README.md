# Pan-cancer analysis of longitudinal metastatic tumors reveals genomic alterations and immune landscape dynamics associated with pembrolizumab sensitivity

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
Serial circulating tumor DNA (ctDNA) monitoring is emerging as a non-invasive strategy to predict and monitor immune checkpoint blockade (ICB) therapeutic efficacy across cancer types. Yet, limited data exist to show the relationship between ctDNA dynamics and tumor genome and immune microenvironment in patients receiving ICB. Here, we present an in-depth analysis of clinical, whole-exome, transcriptome, and ctDNA profiles of 73 patients with advanced solid tumors, across 30 cancer types, from a phase II basket clinical trial of pembrolizumab (NCT02644369). Patients stratified by ctDNA and tumor burden dynamics correspond with survival and clinical benefit. High mutation burden, high expression of immune signatures, and mutations in BRCA2 are associated with pembrolizumab molecular sensitivity, while abundant copy-number alterations and B2M loss-of-heterozygosity corresponded with resistance. Upon treatment, induction of genes expressed by T cell, B cell, and myeloid cell populations are consistent with sensitivity and resistance. We identified the upregulated expression of PLA2G2D, an immune-regulating phospholipase, as a potential biomarker of adaptive resistance to ICB. Together, these findings provide insights into the diversity of immunogenomic mechanisms that underpin pembrolizumab outcomes.

# Repo Contents
This repo contains custom R scripts to reproduce figures presented Yang et al. Nature Communications (2021). Analyses include: 
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
