#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=drift_Fst
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/02/22                  Last Modified: 01/07/24 ###
###########################################################################
###########################################################################
###                     03_drift_fst.sh            						###
###########################################################################

cd $SLURM_SUBMIT_DIR

# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/geneticLoad"
REFDIR="/fs/ess/scratch/PAS1533/smathur/genRes/ref/"

# Get sites with MAF > 0.75

# cd $OUTDIR
# for pop1 in BERG BPNP CCRO CEBO GRLL JENN KBPP KLDR MOSQ PRDF ROME SPVY SSSP WLRD
# do
#     cat freq/$pop1.deleterious.freq | awk '$7 > 0.75 {print $1,$2}' \
#     > drift/$pop1.del.maf.0.75.sites

#     for pop2 in BERG BPNP CCRO CEBO GRLL JENN KBPP KLDR MOSQ PRDF ROME SPVY SSSP WLRD
# 	do
#         if [[ ${pop1} != ${pop2} ]]
# 		then
# 			for r in {1..30}
# 			do
# 			echo "#!/bin/sh -l
# #SBATCH -A PAS1533
# #SBATCH -N 1
# #SBATCH -n 1
# #SBATCH -t 150:00:00
# #SBATCH --job-name=${pop1}.${pop2}.run${r}
# #SBATCH -e %x
# #SBATCH -o %x

# cd $SLURM_SUBMIT_DIR

# module load vcftools
# module load htslib

# cd $OUTDIR

# vcftools --vcf vcf/final152.deleterious.recode.vcf \
# --weir-fst-pop $MAINDIR/upper10/fst/randSamples/${pop1}.shuf8.run${r}.sample.list \
# --weir-fst-pop $MAINDIR/upper10/fst/randSamples/${pop2}.shuf8.run${r}.sample.list \
# --positions drift/$pop1.del.maf.0.75.sites \
# --out drift/persite/${pop1}.drift-${pop2}.run${r}" \
# > $OUTDIR/jobs/${pop1}-${pop2}.run${r}.fst.sh
# 			done
# 		fi
# 	done 
# done 




#### submit jobs ####

cd $OUTDIR

# submit some jobs at one time

#for pops in BERG BPNP
#for pops in CCRO CEBO 
#for pops in GRLL JENN
#for pops in KBPP KLDR 
#for pops in MOSQ PRDF 
#for pops in ROME SPVY 
#for pops in SSSP WLRD

# for pops in ROME SPVY SSSP WLRD
# do
# 	cd $OUTDIR/jobs
# 	for file in ${pops}*.fst.sh
# 	do
# 		cd $OUTDIR/errs
# 		sbatch $OUTDIR/jobs/$file
# 	done
# done

# Get results

# count1=0
# count2=0

mkdir drift/final_results
cd drift/final_results

mkdir meanFst weightFst

for pop1 in BERG BPNP CCRO CEBO GRLL JENN KBPP KLDR MOSQ PRDF ROME SPVY SSSP WLRD
do
	((count1=count1+1))
	count2=0

	for pop2 in BERG BPNP CCRO CEBO GRLL JENN KBPP KLDR MOSQ PRDF ROME SPVY SSSP WLRD
	do
		((count2=count2+1))
		if [[ ${count2} -gt ${count1} ]]
		then
			if [[ ${pop1} != ${pop2} ]]
			then
				for r in {1..30}
				do
					less $OUTDIR/errs/${pop1}.${pop2}.run${r} | grep "mean Fst" | cut -f2 -d ":" | cut -f2 -d " " \
						> meanFst/${pop1}.${pop2}.run${r}.meanFst

					less  $OUTDIR/errs/${pop1}.${pop2}.run${r} | grep "weighted Fst" | cut -f2 -d ":" | cut -f2 -d " " >  \
						weightFst/${pop1}.${pop2}.run${r}.weightFst
				done
			fi
		fi
	done
done
