# Genomic Compatibility Analysis for Assisted Gene Flow
[![DOI](https://img.shields.io/badge/DOI-10.1111%2Fmec.70014-blue)](https://doi.org/10.1111/mec.70014)

This repository contains scripts and data used for the genomic evaluation of assisted gene flow options in endangered Eastern Massasauga rattlesnakes (*Sistrurus catenatus*). The analysis framework assesses genetic compatibility between potential donor and recipient populations for conservation translocation programs.

**Keywords**: conservation genomics, assisted gene flow, genetic compatibility, population genetics, endangered species, rattlesnake, genomic analysis, GATK pipeline

## Table of Contents

- [Citation](#citation)
- [Research Overview](#research-overview)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Interpretation Guide](#interpretation-guide)
- [Related Resources](#related-resources)
- [Support](#support)
- [Species Conservation](#species-conservation)


## Citation

**Mathur, S., & Gibbs, H. L. (2025).** Genomic evaluation of assisted gene flow options in an endangered rattlesnake. *Molecular Ecology*. [https://doi.org/10.1111/mec.70014](https://doi.org/10.1111/mec.70014)


## Research Overview

This study introduces a novel genomic approach to evaluate genetic compatibility for assisted gene flow by analyzing:
- **Deleterious variants**: Mutations potentially harmful to fitness
- **Adaptive variants**: Genetic variants likely under positive selection
- **Local adaptation**: Population-specific beneficial alleles
- **Genetic load**: Burden of harmful mutations

### Key Research Questions
1. What is the genetic impact of introducing individuals from donor populations to recipient populations?
2. How many novel deleterious and adaptive variants would be introduced?
3. What are the effects on masking/unmasking existing deleterious mutations?
4. What is the potential for outbreeding depression through disruption of local adaptation?

## 📁 Repository Structure

```
genomic_compatibility/
├── metadata/                           # Metadata and annotation files
│   ├── Data1.csv                      # Sample metadata
│   ├── SNP_pos/                       # SNP position files by variant type
│   │   ├── final152.Syn.pos.txt      # Synonymous mutations
│   │   ├── final152.deleterious.pos.txt  # Deleterious mutations (LOF + PROVEAN damaging)
│   │   ├── final152.missense.pos.txt     # Missense mutations
│   │   ├── final152.nonSyn.pos.txt       # Non-synonymous mutations
│   │   └── final152.nonsense.pos.txt     # Loss-of-function mutations
│   ├── Scate_HiC.genes.txt           # Reference genome gene annotations
│   ├── final152.sampleList.csv       # Complete sample information
│   ├── upp10_geneCoord.bed.txt       # Coordinates of top 10% DOS genes (adaptive)
│   └── upper10_geneInfo.csv          # Information on adaptive genes
└── scripts/                           # Analysis pipeline scripts
    ├── 01_preprocess/                 # Raw data preprocessing
    ├── 02_filterNstats/              # Quality filtering and statistics
    ├── 03_functionalVariation/       # Functional annotation and analysis
    └── scatenatus_GA_0910.qmd       # Quarto notebook for downstream analysis
```


## Prerequisites

### Software Requirements
- **GATK** (≥4.0): Variant calling and processing
- **BWA**: Read alignment
- **Trimmomatic**: Adapter removal and quality trimming
- **SnpEff**: Functional annotation
- **PROVEAN**: Pathogenicity prediction
- **Samtools**: BAM file manipulation
- **Bcftools**: VCF file processing
- **R** (≥4.0): Statistical analysis and visualization
- **Quarto**: Reproducible document generation

### R Packages
```r
# Required R packages (install before running Quarto notebook)
install.packages(c("tidyverse", "ggplot2", "dplyr", "readr", 
                   "genomics", "vcfR", "adegenet", "hierfstat"))
```


## Interpretation Guide

### Compatibility Metrics

1. **Novel Deleterious Variants**: Number of new harmful mutations introduced to recipient population
2. **Novel Adaptive Variants**: Number of potentially beneficial mutations introduced
3. **Masking Effects**: Number of existing deleterious mutations that become masked (heterozygous)
4. **Unmasking Effects**: Number of existing deleterious mutations that become unmasked (homozygous)

### Decision Framework

**Favorable for assisted gene flow**:
- High introduction of adaptive variants
- Low introduction of deleterious variants
- Positive masking effects
- Minimal disruption of local adaptation

**Unfavorable for assisted gene flow**:
- High genetic load introduction
- Loss of local adaptation
- Excessive outbreeding depression risk

## Related Resources

- **Raw sequence data**: NCBI BioProject [PRJNA1220313](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1220313)
- **Published article**: [Molecular Ecology DOI: 10.1111/mec.70014](https://doi.org/10.1111/mec.70014)
- **GATK Best Practices**: [https://gatk.broadinstitute.org/hc/en-us/sections/360007226651](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)
- **SnpEff Documentation**: [https://pcingola.github.io/SnpEff/](https://pcingola.github.io/SnpEff/)

##  Support

For questions about this analysis framework:
- **Primary contact**: [Samarth Mathur](www.github.com/samarth8392)
- **Lab website**: [Gibbs Lab, Ohio State University](https://u.osu.edu/gibbslab/)
- **Issues**: Please use GitHub Issues for technical problems


## 🐍 Species Conservation

The Eastern Massasauga rattlesnake is listed as **Threatened** under the U.S. Endangered Species Act. This research contributes to evidence-based conservation strategies for this ecologically important species.

---

*⭐ If you find this repository useful for your research, please consider starring it and citing our paper!*

- Conservation genetics community for feedback and collaboration

---

**Keywords**: conservation genomics, assisted gene flow, genetic compatibility, population genetics, endangered species, rattlesnake, genomic analysis, GATK pipeline
