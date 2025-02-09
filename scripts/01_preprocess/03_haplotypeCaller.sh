#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --job-name=haplotypeCaller
#SBATCH -e haplotypeCaller
#SBATCH -o haplotypeCaller

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 07/02/22                  Last Modified: 07/10/22 ###
###########################################################################
###########################################################################
###                     haplotypeCaller.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR

## STEP2: HAPLOTYPECALLER

while read -a sample
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH --partition=gpuserial-40core
#SBATCH --mem=100Gb
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00
#SBATCH --job-name=${sample}_step2
#SBATCH -e %x_%j
#SBATCH -o %x_%j

cd $SLURM_SUBMIT_DIR
module load clara-parabricks
module load picard/2.18.17
module load samtools
module load bwa

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step2/

pbrun haplotypecaller \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam /fs/ess/scratch/PAS1533/smathur/ibd/align/final_bam/${sample}_mark_dups.bam \
--in-recal-file /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/${sample}/${sample}_recalFile.txt \
--gvcf --max-alternate-alleles 2 \
--out-variants ${sample}.raw_variants.g.vcf.gz " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.HapCall.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew2.txt


#### submit jobs ####

while read -a line
do
	cd /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_ind/${line}
	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${line}/${line}.HapCall.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew2.txt


