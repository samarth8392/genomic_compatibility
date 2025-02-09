#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=adapter_removal2
#SBATCH -e adapter_removal2
#SBATCH -o adapter_removal2

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 06/13/22                  Last Modified: 07/08/22 ###
###########################################################################
###########################################################################
###                     adapter_removal.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR

# New sequences 06/13/2022 (85 new + 11 old samples)
# New sequences round 2 07/07/2022 (12 samples)

#Step0: Create folders for jobcodes and errors

#while read -a line
#do
#	mkdir /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${line[0]}
#	mkdir /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_ind/${line[0]}
#done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew2.txt

#run for each sample

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 150:00:00
#SBATCH --job-name=${line[0]}.adptrem
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
#module load fastqc/0.11.8
module load trimmomatic

# Step1: Adapter removal 

java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar \
PE /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/raw/${line[0]}*R1*fastq /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/raw/${line[0]}*R2*fastq \
/fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/${line[0]}.R1.paired /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/${line[0]}.R1.unpaired \
/fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/${line[0]}.R2.paired /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/${line[0]}.R2.unpaired \
LEADING:20 TRAILING:20 MINLEN:30 -threads 20 \
ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:40:10" \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${line[0]}/${line[0]}.adptrem.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew2.txt


while read -a line
do
	cd /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_ind/${line[0]}
	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_ind/${line[0]}/${line[0]}.adptrem.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew2.txt


#cd /fs/ess/scratch/PAS1533/smathur/ibd/data/new22/trimmed/

while read -a line
do
	mv ${line}.R1.paired paired/${line}.R1.filtered.fastq
	mv ${line}.R2.paired paired/${line}.R2.filtered.fastq
done < /fs/ess/scratch/PAS1533/smathur/ibd/align/lists/onlyNew2.txt

# summarize results

#touch adapter_removal.summary

#while read -a line
#do
#	grep "Input Read Pairs:" /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_ind/${line[0]}/${line[0]}.adptrem >> adapter_removal.summary
#done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/newSeq_June22_96.txt

