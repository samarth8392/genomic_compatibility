# Genomic compatibility 🐍
> This GitHub repository contains a collection of scripts developed for the study describing the evaluation of genomic compatibility for assisted gene flow options in an endangered rattlesnake, Eastern Massasauga (_Sistrurus catenatus_).

## Publication
Mathur, S., & Gibbs, H. L. (2025).Genomic evaluation of assisted gene flow options in an endangered rattlesnake. *Molecular Ecology*. *Accepted*. https://doi.org/10.1111/mec.70014



This repository contains the following directories:
> Enter each repo to see more details.

```bash

genomic_compatibility
├── metadata # All the metadata associated with the analysis
│   ├── Data1.csv
│   ├── SNP_pos # SNP positions for various annotated variants
│   │   ├── final152.Syn.pos.txt # Synonymous mutations
│   │   ├── final152.deleterious.pos.txt # Deleterious mutations (LOF + PROVEAN damaging)
│   │   ├── final152.missense.pos.txt # Missense mutations
│   │   ├── final152.nonSyn.pos.txt # Non-Synonymous mutations
│   │   └── final152.nonsense.pos.txt # LoF mutations
│   ├── Scate_HiC.genes.txt # Reference genome genes
│   ├── final152.sampleList.csv # Sample List
│   ├── upp10_geneCoord.bed.txt # Genes in the upper 10% of DOS distribution (Adaptive genes)
│   └── upper10_geneInfo.csv #  Adpative genes info
└── scripts
    ├── 01_preprocess
    │   ├── 01_adapter_removal.sh
    │   ├── 02_alignment.sh
    │   ├── 03_haplotypeCaller.sh
    │   ├── 04_genotypeGVCF.sh
    │   └── 05_gatherVCF.sh
    ├── 02_filterNstats
    │   ├── 01_mapstats.sh
    │   ├── 02_filterVCF.sh
    │   └── 03_vcfStats.sh
    ├── 03_functionalVariation
    │   ├── 01_snpEff.sh
    │   ├── 02_provean.sh
    │   ├── 03_genoFreq.sh
    │   ├── 04_upper10_pi.sh
    │   └── 05_geneDesert.sh
    └── scatenatus_GA_0910.qmd # Quarto cookbook for downstream analysis and data viz


```