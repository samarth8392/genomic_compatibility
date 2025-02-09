#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=phasing
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 11/26/22                 Last Modified: 12/25/22  ###
###########################################################################
###########################################################################
###                     phasing.sh                        				###
###########################################################################


cd $SLURM_SUBMIT_DIR

module load vcftools
module load samtools
module load picard/2.3.0
module load htslib
module load bcftools

### (A) Using SHAPEIT v2 ###

# Copy BAMfiles

#cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/data/

#while read -a sample
#do
#	cp /fs/ess/scratch/PAS1533/smathur/ibdPreNov22/align/scatMap/final_bam/${sample}*bam ./
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/onlyScat_noIOWA.samples.list

#Breaking VCF by chromosome

#while read -a chr
#do
#	vcftools --gzvcf /fs/ess/scratch/PAS1533/smathur/ibd/vcf/all193.allchr.finalSNPs.vcf.gz \
#	--chr ${chr} \
#	--remove /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/STER.final202.sample.list \
#	--recode --recode-INFO-all \
#	--out onlyScat.${chr}
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt


# Creating BAMfiles for each chromosome

#The text file contains three columns:
#The first column gives the sample ID corresponding to the BAM file. 
#The second gives the absolute path of the BAM file. 
#The third gives the "Reference sequence" (chromosome name)

#while read -a chr
#do
#	touch ${chr}.bamlist
#	while read -a sample
#	do
#		echo "${sample} ${sample}_mark_dups.bam ${chr}" >> ${chr}.bamlist
#	done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/onlyScat_noIOWA.samples.list
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt

# Extracting information from sequence reads

while read -a chr
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 150:00:00
#SBATCH --job-name=${chr}.gone
#SBATCH -e %x
#SBATCH -o %x


cd $SLURM_SUBMIT_DIR
module load vcftools

cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/data/

~/softwares/shapeit.v2.904/extractPIRs.v1.r68/extractPIRs \
--bam ${chr}.bamlist \
--vcf onlyScat.${chr}.recode.vcf \
--base-quality 20 --read-quality 20 \
--out onlyScat.${chr}.PIRsList" \
> /fs/ess/scratch/PAS1533/smathur/ibd/phasing/jobs/${chr}.extPIR.sh

done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt

# submit jobs

#while read -a chr
#do
#	cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/errors/
#	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/phasing/jobs/${chr}.extPIR.sh
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt


# Running shape it using read information to phase genotypes into haplotypes

cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/data/

for run in {1..5}
do
	while read -a chr
	do
		echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 150:00:00
#SBATCH --job-name=${chr}_run${run}_phase
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR


cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/data/

~/softwares/shapeit.v2.904/bin/shapeit -assemble --thread 30 \
--input-vcf onlyScat.${chr}.recode.vcf \
--input-pir onlyScat.${chr}.PIRsList \
--states 400 --burn 20 --prune 20 --main 100 --force \
-O onlyScat.phase.${chr}_run${run}" > /fs/ess/scratch/PAS1533/smathur/ibd/phasing/jobs/${chr}_${run}.phase.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt
done


#for run in {1..5}
#do
#	while read -a chr
#	do
#		cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/errors/
#		sbatch /fs/ess/scratch/PAS1533/smathur/ibd/phasing/jobs/${chr}_${run}.phase.sh
#	done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt
#done

# Converting Haplotypes into phased VCF files

cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/haps

#for run in {1..5}
#do
#	while read -a chr
#	do
#		~/softwares/shapeit.v2.904/bin/shapeit -convert \
#        --input-haps onlyScat.phase.${chr}_run${run} \
#        --output-vcf onlyScat.phase.${chr}_run${run}.phased.vcf
#    done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt
#done


# Combine VCFs

# step1: get VCF metadata


#cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/vcf
#
#grep "#" /fs/ess/scratch/PAS1533/smathur/ibd/vcf/all193.allchr.final.SNPs.vcf > vcf.header.txt
#
#for run in {1..5}
#do
#	while read -a chr
#	do
#		cat vcf.header.txt onlyScat.phase.${chr}_run${run}.phased.vcf > vcfWcontigs/onlyScat.wHeader.${chr}_run${run}.phased.vcf
#	done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt
#done


#for run in {1..5}
for run in 1 4
do
	cd /fs/ess/scratch/PAS1533/smathur/ibd/phasing/vcf/vcfWcontigs

	java -jar /usr/local/picard/picard-tools-2.3.0/picard.jar GatherVcfs \
	I=onlyScat.wHeader.Scate-ma1_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-ma2_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-ma3_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-Z_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-ma4_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-ma5_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-ma6_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-ma7_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi1_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi2_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi3_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi4_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi5_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi6_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi7_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi8_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi9_run${run}.phased.vcf \
	I=onlyScat.wHeader.Scate-mi10_run${run}.phased.vcf \
	O=/fs/ess/scratch/PAS1533/smathur/ibd/vcf/phased/onlyScat.allchr.phasedSNPs.run${run}.vcf

	cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf/phased/
	tabix -p vcf onlyScat.allchr.phasedSNPs.run${run}.vcf
done

