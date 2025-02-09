#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=genotypeGVCF
#SBATCH -e %x_%j
#SBATCH -o %x_%j

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 07/10/22                  Last Modified: 08/22/22 ###
###########################################################################
###########################################################################
###                     genotypeGVCF.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module load gatk #gatk/4.1.2.0
#module load vcftools
#module load clara-parabricks


## STEP3: genotypeGVCF

# Using GenomicsDBImport option
# Step1: Create dabtabase of ALL samples for each chromosome separately
# Step2: Do joint genotype calling on ALL samples for each chromosomes
# Step2: Combine vcf of each chr into a single vcf  


# Step1

cd /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/

while read -a chr
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 150:00:00
#SBATCH --job-name=gDB_${chr}
#SBATCH -e %x_%j
#SBATCH -o %x_%j

cd $SLURM_SUBMIT_DIR
module load gatk

cd /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/
mkdir tmp_${chr}
gatk --java-options \"-Xmx170g -XX:+UseParallelGC -XX:ParallelGCThreads=48\" GenomicsDBImport \
--genomicsdb-workspace-path all202.${chr}.gDB \
--batch-size 50 \
-L ${chr} \
--sample-name-map all202.sampleMap \
--tmp-dir=tmp_${chr} \
--reader-threads 10

gatk --java-options \"-Xmx170g -XX:+UseParallelGC -XX:ParallelGCThreads=48\" GenotypeGVCFs \
-R /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
-V gendb://all202.${chr}.gDB \
-O all202.${chr}.vcf 

rm tmp_${chr}" \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_chr/${chr}.gDB.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/chrList.txt
  

#### submit jobs ####

while read -a chr
do
	cd /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_chr/
	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_chr/${chr}.gDB.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/gatk/chrList.txt
