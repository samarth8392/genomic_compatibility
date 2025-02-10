#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=filterVCF
#SBATCH -e filterVCF
#SBATCH -o filterVCF

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 08/25/22                  Last Modified: 11/25/22 ###
###########################################################################
###########################################################################
###                     filterVCF.sh                      				###
###########################################################################

cd $SLURM_SUBMIT_DIR
module load htslib
module load samtools
module load bcftools
module load vcftools

MAINDIR="/fs/ess/scratch/PAS1533/smathur/ibd"

# Filter VCF for each chromosome independently

while read -a chr
do
echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 150:00:00
#SBATCH --job-name=filter_${chr}
#SBATCH -e %x_%j
#SBATCH -o %x_%j

cd $SLURM_SUBMIT_DIR
module load gatk
module load vcftools

cd $MAINDIR/variantCall/gatk/raw/

# Select SNPs

#gatk --java-options \"-Xmx50g -XX:+UseParallelGC -XX:ParallelGCThreads=10\" SelectVariants \
#-R $MAINDIR/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V final152.${chr}.vcf \
#--select-type-to-include SNP \
#-O final152.${chr}.rawSNPs.vcf


#Filter SNPs
# From grossen et al 2020 QD < 2.0, FS > 40.0, SOR > 5.0, MQ< 20.0, −3.0 > MQRandkSum >3.0, 
# −3.0 > ReadPosRankSum >3.0 and AN < 62 (80% of all Alpine ibex individuals)

# using AN < 323 (80% of all samples (2X152 alleles))

#gatk --java-options \"-Xmx50g\" VariantFiltration \
#-R $MAINDIR/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V final152.${chr}.rawSNPs.vcf \
#-filter \"QD < 2.0\" --filter-name \"QD2\" \
#-filter \"MQ < 20.0\" --filter-name \"MQ20\" \
#-filter \"MQRankSum < -3.0\" --filter-name \"MQRankSum-3.0\" \
#-filter \"MQRankSum > 3.0\" --filter-name \"MQRankSum3.0\" \
#-filter \"ReadPosRankSum < -3.0\" --filter-name \"ReadPosRankSum-3.0\" \
#-filter \"ReadPosRankSum > 3.0\" --filter-name \"ReadPosRankSum3.0\" \
#-filter \"SOR > 5.0\" --filter-name \"SOR5\" \
#-filter \"FS > 40.0\" --filter-name \"FS40\" \
#-filter \"AN < 323.0\" --filter-name \"AN323\" \
#-O final152.${chr}.fillterflag_SNPs.vcf

#gatk SelectVariants \
#-R $MAINDIR/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V final152.${chr}.fillterflag_SNPs.vcf \
#-select 'vc.isNotFiltered()' \
#-O final152.${chr}.filltered.SNPs.vcf

# keep only biallelic sites

#vcftools --gzvcf final152.${chr}.filltered.SNPs.vcf \
#--recode --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 \
#--out ../final/final152.${chr}.filtered.final.SNPs.vcf" \
> $MAINDIR/jobcodes/per_chr/${chr}.filterVCF.sh
done < $MAINDIR/variantCall/gatk/chrList2.txt # <- autosomes + unPlaced

#### submit jobs ####

while read -a chr
do
	cd $MAINDIR/errors/per_chr/
	sbatch  $MAINDIR/jobcodes/per_chr/${chr}.filterVCF.sh
done < $MAINDIR/variantCall/gatk/chrList2.txt


### Then follow gatherVCF.sh to combine each chr vcf into single VCF file

### Remove singletons and private doubletons

#get sites
#cd $MAINDIR/sites
vcftools --vcf $MAINDIR/variantCall/gatk/final/final152.allchr.filtered.SNPs.vcf \
--singletons \
--out final152.allchr

cut -f1,2 final152.allchr.singletons > singletons.sites

cd $MAINDIR/vcf
vcftools --vcf $MAINDIR/variantCall/gatk/final/final152.allchr.filtered.SNPs.vcf \
--recode --recode-INFO-all \
--exclude-positions $MAINDIR/sites/singletons.sites \
--out final152.allchr.final.SNPs 

