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
###     Date Created: 06/17/22                  Last Modified: 07/08/22 ###
###########################################################################
###########################################################################
###                     aligment.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR


###############################################################################

##### Sample information ####

# Samples were sequenced in 3 rounds: 2020 (old20), 2021 (old21), 2022 (new22)
# In first round old20, multiple samples were sequenced in different batches (1N,2N,1P,2P,3P,4P,HA)


# List1: Samples with only old20 seq (N=34)
# a) 1batch: 1P/2P/4P/HA ( n = 33 )
samples1_a=( sca0678 sca0689 sca0700 sca0708 sca0742 sca0744 
	sca0807 sca0808 sca0809 sca0884 sca0887 sca0925 sca0926 sca0931 sca0933 sca0982
	sca1037 sca1453 sca1529 sca1530 sca1535 sca1558 sca1559 ste0111 sca1038 sca1057
	ste0121 sca0979 sca1460 ste0005)

# b) 2batch: 2P + HA ( n = 1 )
samples1_b=( sca0886 )

# List2: Samples with both old20 and old21 seq (N=80)
# a) 1batch: 1P/2P (n = 61)
samples2_a=( sca0693 sca0694 sca0709 sca0729 sca0750 sca0802 sca0810 sca0888 
	sca0924 sca0968 sca1027 sca1030 sca1078 sca1541 ste0096 ste0097 ste0098 
	ste0105 ste0120 sca0732 sca0885 sca0965 sca1023 sca1029 sca1036 sca1055 sca1058 sca1272
	sca1372 sca1395 sca1531 sca1536 sca1543 ste0114 sca0673 sca0461 sca0465 sca0472
	sca0529 sca0549 sca0600 sca0604 sca0636 sca0651 sca0901 sca0903 sca0905 sca0917 sca0919 
	sca1223 sca1224 sca1227 sca1228 sca1229 sca1230 sca1231 sca1266 sca1515 )

# b) 2batch: 1N + 2N (n = 17)
samples2_b=( sca0142 sca0143 sca0186 sca0226 sca0236 sca0254 sca0260 sca0863 sca0893
	sca0939 sca0940 sca1018 sca1019 sca1020 sca1026 sca1035 sca1050 )

# c) 2batch: 1P + 3P (n = 2)
samples2_c=( sca0724 sca0812 sca0978 sca1028 )

# d) 2batch: 2P + 4P (n = 2)
samples2_d=( sca1381 )

# List3: Samples with both old20 and new22 seq (N=4)
# a) 1 batch 1P (n = 4)
samples3_a=( sca0734 sca0801 sca0967 sca1542 sca0707 sca0966 sca1532)

# List4: Samples with old20 old21 and new22 seq (N=12)
# a) 1 batch 1P (n = 9)
samples4_a=( sca0531 sca0713 sca0932 sca0962 sca0980 sca1031 sca1259 sca1379 ste0119 )

# c) 2batch: 1N + 2N (n = 3)
samples4_b=( sca0144 sca0180 sca0242 )

#List5: Samples with only new22 seq (N=80)

#################################################################################

#module load bwa
#module load samtools/1.8
#module load picard/2.18.17
#module load vcftools
# step0: Index reference and dictionary


#cd /fs/ess/scratch/PAS1533/smathur/parabricks/data/ref/

#samtools faidx cho_scaff_3genom.fasta
#bwa index cho_scaff_3genom.fasta

#java -jar /apps/picard/2.18.17/picard.jar CreateSequenceDictionary \
#reference=cho_scaff_3genom.fasta \
#output=cho_scaff_3genom.dict


# The Parabricks pipeline:
# LINK: https://docs.nvidia.com/clara/parabricks/v3.5/text/software_overview.html

# (A) Alignment using BWA mem
# (B) Preprocessing: Merge BAM, Co-ordinate sorting, Picard MarkDups, BQSR
# (C) Variant Calling using Haplotype Caller/Deepvariant (NOTE: Try multiple runs with/without different Callers and use union)
# (D) Joint genotyping: GenotypeGVCF
# (E) Variant filtration: Select Variations and Variant Filteration
# (F) Quality checking (QCs)

# Step0: Make folders

#cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/

#while read -a line
#do 
#	mkdir ${line}
#	mkdir ${line}/temp
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/all207_samples.txt

# Get known sites
# Previous VCF thinned for AF > 0.3 and AF < 0.9 ( N = 2,825,863 )

#vcftools --vcf /fs/ess/scratch/PAS1533/smathur/wgr/vcf/onlyScat.10x.wholegenome.SNPs.recode.vcf \
#--maf 0.3 --max-maf 0.9 \
#--recode --recode-INFO-all \
#--out /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites

#################################################################################


# Make scripts for each different list


### samples1_a : only old20 seq  (1 batch) ###

for sample in "${samples1_a[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*fastq \
--out-bam ${sample}_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done 


### samples1_b : only old20 seq  (2 batch 2P/HA) ###

for sample in "${samples1_b[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*2P*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*2P*fastq \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*HA*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*HA*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done 


### samples2_a : both old20 and old21 seq (1 batch) ###

for sample in "${samples2_a[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*fastq \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R2.*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done

### samples2_b : both old20 and old21 seq (2 batch 1N/2N) ###

for sample in "${samples2_b[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*1N*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*1N*fastq \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*2N*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*2N*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R2.*fastq \
--out-bam ${sample}_lib3_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib3

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
I=${sample}_lib3_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done

### samples2_c : both old20 and old21 seq (2 batch 1P/3P) ###

for sample in "${samples2_c[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*1P*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*1P*fastq \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*3P*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*3P*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R2.*fastq \
--out-bam ${sample}_lib3_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib3

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
I=${sample}_lib3_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done

### samples2_d : both old20 and old21 seq (2 batch 2P/4P) ###

for sample in "${samples2_d[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*2P*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*2P*fastq \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*4P*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*4P*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R2.*fastq \
--out-bam ${sample}_lib3_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib3

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
I=${sample}_lib3_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done

### samples3_a : both old20 and new22 seq (1 batch) ###

for sample in "${samples3_a[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*fastq \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R2.*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done

### samples4_a : all 3 (1 batch) ###

for sample in "${samples4_a[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*fastq \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R2.*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R2.*fastq \
--out-bam ${sample}_lib3_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib3

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
I=${sample}_lib3_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done

### samples4_b : all 3 (2 batch 1N/2N) ###

for sample in "${samples4_b[@]}"
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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*1N*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*1N*fastq  \
--out-bam ${sample}_lib1_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R1.*2N*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old20/${sample}.R2.*2N*fastq \
--out-bam ${sample}_lib2_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib2

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/old21/${sample}.R2.*fastq \
--out-bam ${sample}_lib3_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib3 

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R2.*fastq \
--out-bam ${sample}_lib4_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib4

java -jar /apps/picard/2.18.17/picard.jar MergeSamFiles \
I=${sample}_lib1_mark_dups.bam \
I=${sample}_lib2_mark_dups.bam  \
I=${sample}_lib3_mark_dups.bam  \
I=${sample}_lib4_mark_dups.bam  \
O=${sample}_mark_dups.bam


pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done


### samples5 : only new seq ###

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

cd /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/
cd ${sample}

pbrun fq2bam \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-fq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R1.*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/paired/${sample}.R2.*fastq \
--out-bam ${sample}_mark_dups.bam \
--tmp-dir temp \
--read-group-sm ${sample} \
--read-group-lb lib1 

pbrun bqsr \
--ref /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--in-bam ${sample}_mark_dups.bam \
--knownSites /fs/ess/scratch/PAS1533/smathur/ibd/ref/old_vcf/onlyScat.10x.wgr.knownSites.recode.vcf \
--out-recal-file ${sample}_recalFile.txt " \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${sample}/${sample}.align.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew.txt


#### submit jobs ####

while read -a line
do
	cd /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_ind/${line}
	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${line}/${line}.align.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew2.txt


#cd /fs/ess/scratch/PAS1533/smathur/ibd/align/

#touch align_summary3.txt

#while read -a line
#do
#	echo -e "${line}\t" >> align_summary3.txt
#	grep "error" /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_ind/${line}/*_step1_1187* >> align_summary3.txt
#	echo -e "\n" >>  align_summary3.txt
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/all210_samples.txt



#while read -a line
#do
#	mv /fs/ess/scratch/PAS1533/smathur/ibd/align/step1/${line}/${line}_mark_dups.bam /fs/ess/scratch/PAS1533/smathur/ibd/align/final_bam/
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/all210_samples.txt
