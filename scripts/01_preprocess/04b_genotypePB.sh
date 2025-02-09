#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH --partition=gpuserial-40core
#SBATCH --mem=20Gb
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00
#SBATCH --job-name=pb_genomicsDB
#SBATCH -e %x_%j
#SBATCH -o %x_%j

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 07/10/22                  Last Modified: 08/22/22 ###
###########################################################################
###########################################################################
###                     genotypePB.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module load clara-parabricks
#module load picard/2.18.17
#module load samtools
#module load bwa

## STEP3: GENOTYPE GVCF
# Each individual separately

while read -a sample
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH --partition=gpuserial-40core
#SBATCH --mem=20Gb
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00
#SBATCH --job-name=${sample}_genotype
#SBATCH -e %x_%j
#SBATCH -o %x_%j

cd $SLURM_SUBMIT_DIR
module load clara-parabricks
#module load picard/2.18.17
#module load samtools
#module load bwa

cd /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/pbrun/


pbrun genotypegvcf \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-gvcf /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gvcf/${sample}.raw_variants.g.vcf.gz \
--out-vcf ${sample}.raw_variants.vcf " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.genoCall.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/final202.samples.list

#### submit jobs ####

#while read -a line
#do
#	cd /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_ind/${line}
#	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${line}/${line}.genoCall.sh
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/final202.samples.list


### By creating a genomics DB
# see: https://docs.nvidia.com/clara/parabricks/v3.0/text/joint_calling.html

# create genomicsDB
cd /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/pbrun2/

pbrun creategenomicsdb --dir final202.pbrun.genomicsDB

#while read -a sample
#do
#	pbrun importgvcftodb --db-dir final202.pbrun.genomicsDB \
#	--num-threads 5 \
#	--in-gvcf /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gvcf/${sample}.raw_variants.g.vcf.gz
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/final202.samples.list

