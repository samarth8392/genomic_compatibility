#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --job-name=mapStats
#SBATCH -e mapStats
#SBATCH -o mapStats

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 06/17/22                  Last Modified: 07/09/22 ###
###########################################################################
###########################################################################
###                     aligment.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR
MAINDIR="/fs/ess/scratch/PAS1533/smathur/ibd"

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH --job-name=${line}_mapStats
#SBATCH -e ${line}_mapStats
#SBATCH -o ${line}_mapStats

cd $SLURM_SUBMIT_DIR
module load samtools

cd $MAINDIR/align/stats/

samtools depth -a $MAINDIR/align/final_bam/${line}_mark_dups.bam \
| awk '{c++;s+=\$3}END{print s/c}' \
> ${line}.meandepth.txt

samtools depth -a $MAINDIR/align/final_bam/${line}_mark_dups.bam \
| awk '{c++; if(\$3>0) total+=1}END{print (total/c)*100}' \
> ${line}.1xbreadth.txt

samtools depth -a $MAINDIR/align/final_bam/${line}_mark_dups.bam \
| awk '{c++; if(\$3>=10) total+=1}END{print (total/c)*100}' \
> ${line}.10xbreadth.txt

samtools flagstat $MAINDIR/align/final_bam/${line}_mark_dups.bam \
> ${line}.mapstats.txt" \
> $MAINDIR/jobcodes/per_ind/${line}/${line}.mapStats.sh
done < $MAINDIR/lists/final152.sampleList.txt

#### submit jobs ####

while read -a line
do
	cd $MAINDIR/errors/per_ind/${line}
	sbatch  $MAINDIR/jobcodes/per_ind/${line}/${line}.mapStats.sh
done <  $MAINDIR/lists/final152.sampleList.txt
