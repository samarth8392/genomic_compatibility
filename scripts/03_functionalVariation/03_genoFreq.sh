#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=geneticLoad
#SBATCH -e %x
#SBATCH -o %x

###########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     03_genoFreq.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load vcftools


# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes"
SNPEFFDIR="/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/snpEff"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/geneticLoad"

# Get genotypes and frequencies of deleterious and Syn mutations

cd $OUTDIR
mkdir freq geno


for type in deleterious Syn nonsense missense
do
	while read -a pop
	do
		vcftools --vcf $MAINDIR/vcf/onlyScat.final152.phasedSNPs.vcf \
		--positions $MAINDIR/sites/final152.$type.pos.txt \
		--freq \
		--keep $MAINDIR/lists/byPop/${pop}.final152.sample.list \
		--out freq/${pop}.${type}

		sed -i '1d' freq/${pop}.${type}.frq
		less freq/${pop}.${type}.frq | cut -f6 | cut -f2 -d ":" > freq/${pop}.${type}.DAF
		paste freq/${pop}.${type}.frq freq/${pop}.${type}.DAF > freq/${pop}.${type}.freq

		rm freq/*.frq
		rm freq/*.DAF

		vcftools --vcf $MAINDIR/vcf/onlyScat.final152.phasedSNPs.vcf \
		--positions $MAINDIR/sites/final152.$type.pos.txt \
		--extract-FORMAT-info GT \
		--keep $MAINDIR/lists/byPop/${pop}.final152.sample.list \
		--out geno/${pop}.${type}
	done < $MAINDIR/lists/final152.onlyScat.popNames.txt
done


# Sampling N = 8 from each population with sample size > 8 and recalculating derived allele frequency (DAF)

cd /fs/ess/scratch/PAS1533/smathur/ibd/geneticLoad/sample

for pop in BPNP CCRO GRLL JENN KLDR MOSQ PRDF ROME SPVY SSSP WLRD
do
	#mkdir ${pop}
	cd ${pop}
	for type in LoF missense syn del.LoF nondeleterious
	do
		for r in {1..100}
		do
			shuf -n 8 /fs/ess/scratch/PAS1533/smathur/ibd/lists/byPop/${pop}.final152.sample.list > samples/${pop}.shuf8.run${r}.sample.list

			vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/geneticLoad/vcf/onlyScat.final152.${type}.phased.recode.vcf \
			--freq \
			--keep samples/${pop}.shuf8.run${r}.sample.list \
			--out freq/${pop}.${type}.run${r}

			sed -i '1d' freq/${pop}.${type}.run${r}.frq
			less freq/${pop}.${type}.run${r}.frq | cut -f6 | cut -f2 -d ":" > freq/${pop}.${type}.run${r}.DAF
			paste freq/${pop}.${type}.run${r}.frq freq/${pop}.${type}.run${r}.DAF > freq/${pop}.${type}.run${r}.freq

			rm freq/*.frq
			rm freq/*.DAF

			vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/geneticLoad/vcf/onlyScat.final152.${type}.phased.recode.vcf \
			--extract-FORMAT-info GT \
			--keep samples/${pop}.shuf8.run${r}.sample.list \
			--out geno/${pop}.${type}.run${r}

		done
	done
	cd /fs/ess/scratch/PAS1533/smathur/ibd/geneticLoad/sample
done



