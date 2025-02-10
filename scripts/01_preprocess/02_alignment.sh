#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --job-name=aligment4
#SBATCH -e aligment4
#SBATCH -o aligment4

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     aligment.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR

MAINDIR="/fs/ess/scratch/PAS1533/smathur/ibd"


while read -a sample
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH --partition=gpuserial-40core
#SBATCH --mem=100Gb
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00
#SBATCH --job-name=${sample}_step1
#SBATCH -e %x_%j
#SBATCH -o %x_%j


cd $SLURM_SUBMIT_DIR
module load clara-parabricks
module load picard/2.18.17
module load samtools
module load bwa


## STEP1: ALIGNMENT AND PREPROCESSING

cd $MAINDIR/align/step1/
cd ${sample}

pbrun fq2bam \
--ref $MAINDIR/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq $MAINDIR/data/new22/trimmed/paired/${sample}.R1.*fastq $MAINDIR/data/new22/trimmed/paired/${sample}.R2.*fastq \
--out-bam ${sample}_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun bqsr \
--ref $MAINDIR/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites $MAINDIR/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> $MAINDIR/jobcodes/per_ind/${sample}/${sample}.align.sh
done < $MAINDIR/lists/final152.sampleList.txt


#### submit jobs ####

while read -a line
do
	cd $MAINDIR/errors/per_ind/${line}
	sbatch  $MAINDIR/jobcodes/per_ind/${line}/${line}.align.sh
done < $MAINDIR/align/lists/final152.sampleList.txt
