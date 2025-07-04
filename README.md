# 🧬 Genomic Compatibility for Assisted Gene Flow

> **Evaluating genomic compatibility for conservation of the endangered Eastern Massasauga rattlesnake**

[![DOI](https://img.shields.io/badge/DOI-10.1111%2Fmec.70014-blue)](https://doi.org/10.1111/mec.70014)


## 📋 Table of Contents

- [🧬 Genomic Compatibility for Assisted Gene Flow](#-genomic-compatibility-for-assisted-gene-flow)
  - [📋 Table of Contents](#-table-of-contents)
  - [🐍 About](#-about)
  - [📖 Publication](#-publication)
  - [🔬 Key Features](#-key-features)
  - [🗂️ Repository Structure](#️-repository-structure)
    - [Scripts](#scripts)
  - [🔍 Key Analyses](#-key-analyses)
    - [🧪 Functional Variant Classification](#-functional-variant-classification)
    - [📊 Population Genomics Metrics](#-population-genomics-metrics)
    - [🎯 Adaptive Genomics](#-adaptive-genomics)
  - [📈 Results Highlights](#-results-highlights)
  - [📧 Contact](#-contact)
  - [🐍 Species Conservation](#-species-conservation)

## 🐍 About

This repository contains the complete computational pipeline for assessing genomic compatibility in assisted gene flow strategies for the **Eastern Massasauga rattlesnake** (*Sistrurus catenatus*) - a critically endangered species facing habitat fragmentation and population decline.

Our genomic approach evaluates the genetic consequences of mixing populations to inform evidence-based conservation decisions and optimize genetic rescue efforts.

## 📖 Publication

**Mathur, S., & Gibbs, H. L. (2025).** *Genomic evaluation of assisted gene flow options in an endangered rattlesnake.* **Molecular Ecology.** [https://doi.org/10.1111/mec.70014](https://doi.org/10.1111/mec.70014)

## 🔬 Key Features

- **Comprehensive Variant Analysis**: From raw sequencing data to functional variant annotation
- **Population Genomics**: Multi-population comparison and genetic diversity assessment  
- **Functional Genomics**: Identification of deleterious mutations and adaptive genes
- **Conservation Genomics**: Genomic compatibility metrics for assisted gene flow
- **Reproducible Pipeline**: End-to-end workflow with quality control and visualization

## 🗂️ Repository Structure

```
genomic_compatibility/
├── 📁 metadata/              # Sample information and annotations
│   ├── 📊 Data1.csv          # Primary dataset
│   ├── 📋 final152.sampleList.csv  # Complete sample list
│   ├── 🧬 Scate_HiC.genes.txt      # Reference genome annotations
│   ├── 📁 SNP_pos/           # Variant position files by functional category
│   └── 🎯 upp10_geneCoord.bed.txt  # Adaptive gene coordinates
├── 📁 scripts/               # Analysis pipeline scripts
│   ├── 📁 01_preprocess/     # Data preprocessing & alignment
│   ├── 📁 02_filterNstats/   # Quality control & filtering
│   ├── 📁 03_functionalVariation/  # Functional annotation
│   └── 📊 scatenatus_GA_0910.qmd   # Main analysis & visualization
└── 📄 README.md              # This file
```


### Scripts

1. **Data Preprocessing**
   ```bash
   cd scripts/01_preprocess/
   bash 01_adapter_removal.sh
   bash 02_alignment.sh
   bash 03_haplotypeCaller.sh
   ```

2. **Quality Control & Filtering**
   ```bash
   cd ../02_filterNstats/
   bash 01_mapstats.sh
   bash 02_filterVCF.sh
   ```

3. **Functional Analysis**
   ```bash
   cd ../03_functionalVariation/
   bash 01_snpEff.sh
   bash 02_provean.sh
   ```

4. **Downstream Analysis & Visualization**
   ```bash
   quarto render scatenatus_GA_0910.qmd
   ```

## 🔍 Key Analyses

### 🧪 Functional Variant Classification
- **Synonymous**: Silent mutations (no amino acid change)
- **Missense**: Amino acid substitutions
- **Nonsense**: Premature stop codons
- **Deleterious**: Loss-of-function + PROVEAN damaging variants

### 📊 Population Genomics Metrics
- Nucleotide diversity (π)
- Genetic differentiation (FST)
- Linkage disequilibrium patterns
- Demographic history inference

### 🎯 Adaptive Genomics
- Detection of genes under selection
- Functional enrichment analysis
- Gene desert identification
- Recombination rate variation

## 📈 Results Highlights

- **152 individuals** analyzed across multiple populations
- **Genome-wide variant discovery** with functional annotation
- **Population structure** and genetic diversity characterization
- **Genomic compatibility assessment** for conservation strategies


## 📧 Contact

**Samarth Mathur** 
📧 mathur.112@osu.edu  
🏛️ Ohio State University

**H. Lisle Gibbs** - Principal Investigator  
📧 gibbs.128@osu.edu  
🏛️ Ohio State University

## 🐍 Species Conservation

The Eastern Massasauga rattlesnake is listed as **Threatened** under the U.S. Endangered Species Act. This research contributes to evidence-based conservation strategies for this ecologically important species.

---

*⭐ If you find this repository useful for your research, please consider starring it and citing our paper!*
