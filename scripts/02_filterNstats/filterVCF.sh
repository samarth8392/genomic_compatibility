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

cd /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/raw/

# Select SNPs

#gatk --java-options \"-Xmx50g -XX:+UseParallelGC -XX:ParallelGCThreads=10\" SelectVariants \
#-R /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V all202.${chr}.vcf \
#--select-type-to-include SNP \
#-O all202.${chr}.rawSNPs.vcf


#Filter SNPs
# From grossen et al 2020 QD < 2.0, FS > 40.0, SOR > 5.0, MQ< 20.0, −3.0 > MQRandkSum >3.0, 
# −3.0 > ReadPosRankSum >3.0 and AN < 62 (80% of all Alpine ibex individuals)

# using AN < 323 (80% of all samples (2X202 alleles))

#gatk --java-options \"-Xmx50g\" VariantFiltration \
#-R /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V all202.${chr}.rawSNPs.vcf \
#-filter \"QD < 2.0\" --filter-name \"QD2\" \
#-filter \"MQ < 20.0\" --filter-name \"MQ20\" \
#-filter \"MQRankSum < -3.0\" --filter-name \"MQRankSum-3.0\" \
#-filter \"MQRankSum > 3.0\" --filter-name \"MQRankSum3.0\" \
#-filter \"ReadPosRankSum < -3.0\" --filter-name \"ReadPosRankSum-3.0\" \
#-filter \"ReadPosRankSum > 3.0\" --filter-name \"ReadPosRankSum3.0\" \
#-filter \"SOR > 5.0\" --filter-name \"SOR5\" \
#-filter \"FS > 40.0\" --filter-name \"FS40\" \
#-filter \"AN < 323.0\" --filter-name \"AN323\" \
#-O all202.${chr}.fillterflag_SNPs.vcf

#gatk SelectVariants \
#-R /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V all202.${chr}.fillterflag_SNPs.vcf \
#-select 'vc.isNotFiltered()' \
#-O all202.${chr}.filltered.SNPs.vcf

# keep only biallelic sites

#vcftools --gzvcf all202.${chr}.filltered.SNPs.vcf \
#--recode --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 \
#--out ../final/all202.${chr}.filtered.final.SNPs.vcf" \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_chr/${chr}.filterVCF.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/chrList2.txt

#### submit jobs ####

#while read -a chr
#do
#	cd /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_chr/
#	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_chr/${chr}.filterVCF.sh
#done < /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/chrList2.txt


### Remove singletons and private doubletons

#get sites
#cd /fs/ess/scratch/PAS1533/smathur/ibd/sites
#vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/final/all202.allchr.filtered.SNPs.vcf \
#--singletons \
#--out all202.allchr

#cut -f1,2 all202.allchr.singletons > singletons.sites

#cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf
#vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/final/all202.allchr.filtered.SNPs.vcf \
#--recode --recode-INFO-all \
#--exclude-positions /fs/ess/scratch/PAS1533/smathur/ibd/sites/singletons.sites \
#--out all202.allchr.nosing.SNPs 

# After filtering, kept 32102010 out of a possible 39565443 Sites


### Remove SNPs that are only polymorphic/fixed in Ster (i.e. missing in Scatenatus)

# get sites

cd /fs/ess/scratch/PAS1533/smathur/ibd/sites

#vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/vcf/all202.allchr.nosing.SNPs.vcf \
#--freq \
#--remove /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/STER.final202.sample.list \
#--out onlyScat.SNPs

#cut -f6 onlySter.SNPs.frq | cut -f2 -d ":" > allSNP.sites

#cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf
#vcftools --vcf all202.allchr.nosing.SNPs.vcf \
#--recode --recode-INFO-all \
#--positions /fs/ess/scratch/PAS1533/smathur/ibd/sites/scatSNPs.txt \
#--out all202.allchr.final.SNPs 

# After filtering, kept 19774972 out of a possible 32102010 Sites

#### Remove IOWA samples ####

#cd /fs/ess/scratch/PAS1533/smathur/ibd/sites/

#vcftools --gzvcf /fs/ess/scratch/PAS1533/smathur/ibdPreNov22/vcf/all202.allchr.finalSNPs.vcf.gz \
#--freq \
#--remove /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/IOWA.final202.sample.list \
#--out all202.allchr.noIOWA.SNPs

#cut -f6 all202.allchr.noIOWA.SNPs.frq | cut -f2 -d ":" > all202.allchr.noIOWA.SNPs.freq

cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf
vcftools --gzvcf /fs/ess/scratch/PAS1533/smathur/ibdPreNov22/vcf/all202.allchr.finalSNPs.vcf.gz \
--recode --recode-INFO-all \
--remove /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/IOWA.final202.sample.list \
--positions /fs/ess/scratch/PAS1533/smathur/ibd/sites/all202.allchr.noIOWA.sites \
--out all193.allchr.final.SNPs 

#### After filtering, kept 18995511 out of a possible 19774972 Sites #############

mv all193.allchr.final.SNPs.recode.vcf all193.allchr.final.SNPs.vcf

bgzip -c all193.allchr.final.SNPs.vcf > all193.allchr.finalSNPs.vcf.gz
tabix -p vcf all193.allchr.finalSNPs.vcf.gz
cp all193.allchr.finalSNPs.vcf.gz* /fs/ess/PAS1533/users/smathur/ibd/
