# Genomic compatibility 🐍
> This GitHub repository contains a collection of scripts developed for the study describing the evaluation of genomic compatibility for assisted gene flow options in an endangered rattlesnake, Eastern Massasauga (_Sistrurus catenatus_).

This repository contains the following directories:
> Enter each repo to see more details.

```bash

genomic_compatibility
├── metadata
│   ├── Data1.csv
│   ├── SNP_pos
│   │   ├── final152.Syn.pos.txt
│   │   ├── final152.deleterious.pos.txt
│   │   ├── final152.missense.pos.txt
│   │   ├── final152.nonSyn.pos.txt
│   │   └── final152.nonsense.pos.txt
│   ├── Scate_HiC.genes.txt
│   ├── final152.sampleList.csv
│   ├── upp10_geneCoord.bed.txt
│   └── upper10_geneInfo.csv
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
    └── scatenatus_GA_0910.qmd


```


```