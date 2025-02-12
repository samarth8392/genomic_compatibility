#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=upper10
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     04_upper10.sh              						###
###########################################################################

cd $SLURM_SUBMIT_DIR

# module load vcftools
module load htslib
module load samtools

# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/geneDesert"
REFDIR="/fs/ess/scratch/PAS1533/smathur/genRes/ref/"


####################### 
# GENE DESERT REGIONS #
####################### 


cd $OUTDIR

cd $SLURM_SUBMIT_DIR
module load vcftools
module load bedops

# To generate distances between each SNP and its nearest genes upstream and downstream, you will need a .gff file 
# To prep the .gff file, first remove the trailing semicolons in the file and then remove all features that aren’t 
# genes (e.g., mRNA, exons, etc.)

cd $REFDIR/annotations

awk '{gsub(/;$/,"");print}' Scate_HiC_rnd4.all.putative.function.gff > Scate_HiC.revisedgenome.all.gff
cat Scate_HiC.revisedgenome.all.gff | awk '$3 =="gene" {print $0}' > Scate_HiC.revisedgenome.all_justgenes.gff
cat Scate_HiC.revisedgenome.all.gff | awk '$3 =="CDS" {print $0}' > Scate_HiC.revisedgenome.all_justCDS.gff

# Use bedops to measure distances between genes and SNPs:

cd $OUTDIR

gff2bed < $REFDIR/annotations/Scate_HiC.revisedgenome.all_justgenes.gff > Scate_HiC.revisedgenome.all_justgenes.bed
gff2bed < $REFDIR/annotations/Scate_HiC.revisedgenome.all_justCDS.gff > Scate_HiC.revisedgenome.all_justCDS.bed
vcf2bed < $MAINDIR/vcf/final152.allchr.finalSNPs.vcf > onlyScat.final152.phasedSNPs.bed
sort-bed onlyScat.final152.phasedSNPs.bed > onlyScat.final152.phasedSNPs.sorted.bed
sort-bed Scate_HiC.revisedgenome.all_justgenes.bed > Scate_HiC.sortedgff_justgenes.bed

closest-features --delim '\t' --dist onlyScat.final152.phasedSNPs.sorted.bed Scate_HiC.sortedgff_justgenes.bed \
> onlyScat.final152.phasedSNPs.distance2genes_tab.txt

# Some markers will be upstream from a gene but there will be no genes downstream (or vice versa).  If we ignore this 
# we can end up incorporating a lot of SNPs that are at the extreme edge of a scaffold.  This makes me nervous because 
# for all we know there is a gene 100 bp downstream of the marker, but we can't tell because the scaffold doesn't extend 
# that far.  As such, I prefer to work with SNPs that have an annotated gene both upstream and downstream, so I know 
# exactly what kind of distances I'm working with.

# The following selects for SNPs that have a SNP both upstream and downstream 

cat onlyScat.final152.phasedSNPs.distance2genes_tab.txt | \
    awk '{if ($194!~/NA/) print $0;}' | \
    awk '{if ($195!~/NA/) print $0;}' \
    > distance2genes_noNA_clean.txt

# Keep only useful columns
cat distance2genes_noNA_clean.txt | \
    cut -f1,2,3,194,195,196,203,204,205,206,207,214,215 \
    > distance2genes_clean_all.txt


################################################################################################
# These steps were done in R:

module load R
R
setwd("$OUTDIR") 

df <- read.csv("distance2genes_clean_all.txt", header=F, sep="\t")

# > dim(df)
# [1] 18,603,374       13

# Remove SNPs that are actually within annotated genes.
df.nogene <- df[-which(df$V8 == 0 | df$V13 == 0),] 

# > dim(df.nogene)
# [1] 15,012,745       13

# Remove SNPs that are less than 1Mb on either side of an annotated gene
df.gendes <- df.nogene[which(abs(df.nogene$V8) > 1e6 & abs(df.nogene$V13) > 1e6), ] 

# > dim(df.gendes)
# [1] 175,545     13

snps <- df.gendes[,c(1,2,3)]

write.table(df.gendes, "Snps_in_genedeserts.txt",quote=F,row.names=F, col.names=F)
write.table(snps, "genedeserts_snps.bed",quote=F,row.names=F, col.names=F)
################################################################################################

# Get gene Desert VCFs

cd $OUTDIR

#final 152
vcftools --vcf $MAINDIR/vcf/final152.allchr.finalSNPs.vcf \
--recode --recode-INFO-all \
--keep $MAINDIR/lists/final152.sampleList \
--bed beds/genedeserts_snps.bed --out vcf/onlyScat.final152.geneDes


#onlyOH
vcftools --vcf $MAINDIR/vcf/final152.allchr.finalSNPs.vcf \
--recode --recode-INFO-all \
--keep $MAINDIR/lists/onlyOH.sampleList \
--bed beds/genedeserts_snps.bed --out vcf/onlyScat.onlyOH.geneDes