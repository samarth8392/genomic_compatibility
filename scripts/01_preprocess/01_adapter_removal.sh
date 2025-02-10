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
###########################################################################
###########################################################################
###                     adapter_removal.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR

MAINDIR="/fs/ess/scratch/PAS1533/smathur/ibd"


#Step0: Create folders for jobcodes and errors

#while read -a line
#do
#	mkdir $MAINDIR/jobcodes/per_ind/${line[0]}
#	mkdir $MAINDIR/errors/per_ind/${line[0]}
#done < $MAINDIR/lists/final152.sampleList.txt

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
PE $MAINDIR/data/new22/raw/${line[0]}*R1*fastq $MAINDIR/data/new22/raw/${line[0]}*R2*fastq \
$MAINDIR/data/new22/trimmed/${line[0]}.R1.paired $MAINDIR/data/new22/trimmed/${line[0]}.R1.unpaired \
$MAINDIR/data/new22/trimmed/${line[0]}.R2.paired $MAINDIR/data/new22/trimmed/${line[0]}.R2.unpaired \
LEADING:20 TRAILING:20 MINLEN:30 -threads 20 \
ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:40:10" \
> $MAINDIR/jobcodes/per_ind/${line[0]}/${line[0]}.adptrem.sh
done < $MAINDIR/lists/final152.sampleList.txt


while read -a line
do
	cd $MAINDIR/errors/per_ind/${line[0]}
	sbatch  $MAINDIR/jobcodes/per_ind/${line[0]}/${line[0]}.adptrem.sh
done < $MAINDIR/lists/final152.sampleList.txt


#cd $MAINDIR/data/new22/trimmed/

while read -a line
do
	mv ${line}.R1.paired paired/${line}.R1.filtered.fastq
	mv ${line}.R2.paired paired/${line}.R2.filtered.fastq
done < $MAINDIR/lists/final152.sampleList.txt

# summarize results

#touch adapter_removal.summary

#while read -a line
#do
#	grep "Input Read Pairs:" $MAINDIR/errors/per_ind/${line[0]}/${line[0]}.adptrem >> adapter_removal.summary
#done < $MAINDIR/lists/final152.sampleList.txt

