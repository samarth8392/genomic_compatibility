#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=vcf.final152
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     vcfStats.sh              						###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load vcftools
module load htslib

MAINDIR="/fs/ess/scratch/PAS1533/smathur/ibd"

cd $MAINDIR/vcf/

#### Variant statistics ####


for i in depth site-pi het relatedness2 missing-indv missing-site
do
	vcftools --vcf final152.allchr.final.SNPs.vcf \
	--$i \
	--out stats/final152.allchr
done

# Get vcfs by pop

cd $MAINDIR/vcf/

while read -a pop
do
	vcftools --vcf final152.allchr.final.SNPs.vcf \
	--recode --recode-INFO-all \
	--keep $MAINDIR/lists/bypop/${pop}.final152.sample.list \
	--out bypop/${pop}.final152.allchr.finalSNPs
done < $MAINDIR/lists/popnames.txt

# Get diversity by pop

while read -a pop
do
	vcftools --vcf final152.allchr.final.vcf \
	--window-pi 100000 --window-pi-step 100000 \
	--keep $MAINDIR/lists/bypop/${pop}.final152.sample.list \
	--out $MAINDIR/diversity/bypop/pi/${pop}.100kb

	vcftools --vcf final152.allchr.final.vcf \
	--freq \
	--keep $MAINDIR/lists/bypop/${pop}.final152.sample.list \
	--out $MAINDIR/diversity/bypop/freq/${pop}

	vcftools --vcf final152.allchr.final.vcf \
	--extract-FORMAT-info GT \
	--keep $MAINDIR/lists/bypop/${pop}.final152.sample.list \
	--out $MAINDIR/diversity/bypop/geno/${pop}

done < $MAINDIR/lists/popnames.txt

