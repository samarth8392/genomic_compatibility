#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=vcf.final172
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/02/22                  Last Modified: 01/31/23 ###
###########################################################################
###########################################################################
###                     vcftools.sh              						###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load vcftools
module load htslib

#cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf/

#### Variant statistics ####

### Whole genome (all193 = 18,995,511 )

#for i in depth site-pi het relatedness2 missing-indv missing-site
#do
#	vcftools --vcf all193.allchr.final.SNPs.vcf \
#	--$i \
#	--out stats/all193.allchr
#done

# Get vcfs by pop

#cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf/

#while read -a pop
#do
#	vcftools --vcf all193.allchr.final.SNPs.vcf \
#	--recode --recode-INFO-all \
#	--keep /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/${pop}.final202.sample.list \
#	--out bypop/${pop}.all193.allchr.finalSNPs
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/popnames.txt

# Get diversity by pop

#while read -a pop
#do
	#vcftools --vcf all202.allchr.final.vcf \
	#--window-pi 100000 --window-pi-step 100000 \
	#--keep /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/${pop}.final202.sample.list \
	#--out /fs/ess/scratch/PAS1533/smathur/ibd/diversity/bypop/pi/${pop}.100kb

	#vcftools --vcf all202.allchr.final.vcf \
	#--freq \
	#--keep /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/${pop}.final202.sample.list \
	#--out /fs/ess/scratch/PAS1533/smathur/ibd/diversity/bypop/freq/${pop}

	#vcftools --vcf all202.allchr.final.vcf \
	#--extract-FORMAT-info GT \
	#--keep /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/${pop}.final202.sample.list \
	#--out /fs/ess/scratch/PAS1533/smathur/ibd/diversity/bypop/geno/${pop}

#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/popnames.txt

#cd /fs/ess/scratch/PAS1533/smathur/ibd/pca/

#vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/treemix/all193.auto.noN.finalSNPs.LDpruned.vcf \
#--recode --recode-INFO-all \
#--remove /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/STER.final202.sample.list \
#--out onlyScat.auto.noN.LDpruned

#cd /fs/ess/scratch/PAS1533/smathur/ibd/geneticLoad/phased/vcf

#for type in del.LoF nondeleterious syn
#for type in syn
#do
	#vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/vcf/phased/onlyScat.allchr.phasedSNPs.run1.vcf \
	#--recode --recode-INFO-all \
	#--positions /fs/ess/scratch/PAS1533/smathur/ibd/dnds/all193.${type}.pos.txt \
	#--out all193.${type}.phased

	#while read -a pop
	#do
	#	vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/vcf/phased/onlyScat.allchr.phasedSNPs.run1.vcf \
	#	--recode --recode-INFO-all \
	#	--positions /fs/ess/scratch/PAS1533/smathur/ibd/dnds/all193.${type}.pos.txt \
	#	--keep /fs/ess/scratch/PAS1533/smathur/ibd/lists/bypop/${pop}.final202.sample.list \
	#	--out bypop/${pop}.all193.${type}.phased
	#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/popnames_noIOWA.txt
#done

# del.LoF = After filtering, kept 7770 out of a possible 18686655 Sites
# nondeleterious = After filtering, kept 61035 out of a possible 18686655 Sites
# syn = After filtering, kept 94040 out of a possible 18686655 Sites


################################################################################################


## UPDATE: 01/31/23 (REMOVE GRLL-2/3/4) ##


#cd /fs/ess/scratch/PAS1533/smathur/ibd/sites/

#vcftools --gzvcf /fs/ess/scratch/PAS1533/smathur/ibd/IBDBackupJan23/vcf/all193.allchr.finalSNPs.vcf.gz \
#--freq \
#--keep /fs/ess/scratch/PAS1533/smathur/ibd/lists/final172.samples.list \
#--out final172

#sed -i '1d' final172.frq
#less final172.frq | cut -f6 | cut -f2 -d ":" > final172.DAF
#paste final172.frq final172.DAF > final172.freq
#rm final172.frq
#rm final172.DAF
#cat final172.freq | awk '{if ($7!= 0) print $0;}' | cut -f1,2 > final172.allchr.SNP.sites # N = 18,923,075

#cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf

#vcftools --gzvcf /fs/ess/scratch/PAS1533/smathur/ibd/IBDBackupJan23/vcf/all193.allchr.finalSNPs.vcf.gz \
#--positions /fs/ess/scratch/PAS1533/smathur/ibd/sites/final172.allchr.SNP.sites \
#--keep /fs/ess/scratch/PAS1533/smathur/ibd/lists/final172.samples.list \
#--recode --recode-INFO-all \
#--out final172.allchr.finalSNPs

#mv final172.allchr.finalSNPs.recode.vcf final172.allchr.finalSNPs.vcf

#bgzip -c final172.allchr.finalSNPs.vcf > final172.allchr.finalSNPs.vcf.gz
#tabix -p vcf final172.allchr.finalSNPs.vcf.gz

# After filtering, kept 18,923,075 out of a possible 18995511 Sites



## Phased ##

cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf/phase

for r in {1..5}
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=phase.run${r}
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR

module load vcftools
module load bcftools
module load htslib

cd /fs/ess/scratch/PAS1533/smathur/ibd/vcf/phase

vcftools --vcf /fs/ess/scratch/PAS1533/smathur/ibd/IBDBackupJan23/vcf/phased/onlyScat.allchr.phasedSNPs.run${r}.vcf \
--positions /fs/ess/scratch/PAS1533/smathur/ibd/sites/final172.allchr.SNP.sites \
--recode --recode-INFO-all \
--out onlyScat.final172.phasedSNPs.run${r}

mv onlyScat.final172.phasedSNPs.run${r}.recode.vcf onlyScat.final172.phasedSNPs.run${r}.vcf

tabix -p vcf onlyScat.final172.phasedSNPs.run${r}.vcf" \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/random/phase.final172.run${r}.sh
done

# After filtering, kept 18615434 out of a possible 18686655 Sites
