# Genomic Compatibility Repository File Descriptions

## Scripts

### 01_preprocess Subdirectory

**01_adapter_removal.sh** - Shell script for removing sequencing adapters and low-quality bases from raw FASTQ files.

**02_alignment.sh** - Script for aligning cleaned reads to the reference genome using BWA or similar alignment tools.

**03_haplotypeCaller.sh** - GATK HaplotypeCaller script for variant discovery and genotyping from aligned BAM files.

**04_genotypeGVCF.sh** - Script for joint genotyping of genomic variant call format (GVCF) files across all samples.

**05_gatherVCF.sh** - Utility script for combining and gathering VCF files from multiple chromosomes or regions.

### 02_filterNstats Subdirectory

**01_mapstats.sh** - Script for calculating mapping statistics and quality metrics from aligned BAM files.

**02_filterVCF.sh** - VCF filtering script applying quality thresholds, missing data filters, and population genetic criteria.

**03_vcfStats.sh** - Script for generating comprehensive statistics and summaries from filtered VCF files.

### 03_functionalVariation Subdirectory

**01_snpEff.sh** - Script for functional annotation of variants using SnpEff to predict biological effects of mutations.

**02_provean.sh** - PROVEAN analysis script for predicting the functional impact of amino acid substitutions.

**03_genoFreq.sh** - Script for calculating genotype and allele frequencies across populations and functional variant categories.

**04_upper10_pi.sh** - Analysis script for calculating nucleotide diversity (π) in genes within the upper 10% of the DOS distribution.

**05_geneDesert.sh** - Script for identifying and analyzing gene desert regions with reduced recombination rates.

## Main Analysis File

**scatenatus_GA_0910.qmd** - Quarto markdown document containing the complete downstream analysis pipeline, statistical tests, and data visualizations for the genomic compatibility study.