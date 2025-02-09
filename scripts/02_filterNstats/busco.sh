#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 150:00:00
#SBATCH --job-name=busco
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/02/22                  Last Modified: 11/10/22 ###
###########################################################################
###########################################################################
###                     provean.sh  				                    ###
###########################################################################

cd $SLURM_SUBMIT_DIR

module purge
module load python
conda activate busco

#export AUGUSTUS_CONFIG_PATH=/scratch/bell/mathur20/osu/ibd/busco/augustus/config

cd /fs/ess/scratch/PAS1533/smathur/ibd/busco

busco --in /fs/ess/scratch/PAS1533/smathur/ibd/ref/scat/Scatenatus_HiC_v1.1.fasta \
--augustus -l tetrapoda_odb10/ -m genome -c 16 -f \
-o busco_results 

generate_plot.py -wd busco_results