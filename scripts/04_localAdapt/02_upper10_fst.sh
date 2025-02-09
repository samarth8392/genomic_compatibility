#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=fst
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/02/22                  Last Modified: 01/06/24 ###
###########################################################################
###########################################################################
###                     fst.sh              							###
###########################################################################

###########################################################################

cd $SLURM_SUBMIT_DIR
module load vcftools


# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes/"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/upper10"
VCFDIR="/fs/ess/scratch/PAS1533/smathur/genRes/upper10/vcf"


### Using only N=8 samples from populations with N>8 ###

cd $OUTDIR

# for pop in BPNP CCRO GRLL JENN KLDR MOSQ PRDF ROME SPVY SSSP WLRD
# do
# 	for r in {1..30}
# 	do
# 		shuf -n 8 $MAINDIR/lists/byPop/${pop}.final172.sample.list \
# 		> fst/randSamples/${pop}.shuf8.run${r}.sample.list
# 	done
# done

# for pop in BERG CEBO KBPP
# do
# 	for r in {1..30}
# 	do
# 		cp $MAINDIR/lists/byPop/${pop}.final172.sample.list \
# 		fst/randSamples/${pop}.shuf8.run${r}.sample.list
# 	done
# done


# count1=0
# count2=0

# for pop1 in BERG BPNP CCRO CEBO GRLL JENN KBPP KLDR MOSQ PRDF ROME SPVY SSSP WLRD
# do
# 	((count1=count1+1))
# 	count2=0

# 	for pop2 in BERG BPNP CCRO CEBO GRLL JENN KBPP KLDR MOSQ PRDF ROME SPVY SSSP WLRD
# 	do
# 		((count2=count2+1))
# 		if [[ ${count2} -gt ${count1} ]]
# 		then
# 			if [[ ${pop1} != ${pop2} ]]
# 			then
# 				for r in {1..30}
# 				do
# 					echo "#!/bin/sh -l
# #SBATCH -A PAS1533
# #SBATCH -N 1
# #SBATCH -n 1
# #SBATCH -t 150:00:00
# #SBATCH --job-name=${pop1}.${pop2}.run${r}
# #SBATCH -e %x
# #SBATCH -o %x

# cd $SLURM_SUBMIT_DIR

# module load vcftools

# cd $OUTDIR

# vcftools --vcf $VCFDIR/final152.upper10.nonSyn.phasedSNPs.recode.vcf \
# --weir-fst-pop fst/randSamples/${pop1}.shuf8.run${r}.sample.list \
# --weir-fst-pop fst/randSamples/${pop2}.shuf8.run${r}.sample.list \
# --out fst/${pop1}-${pop2}.shuf8.up10.nonSyn.run${r}" \
# > jobs/${pop1}-${pop2}.shuf8.run${r}.fst.sh
# 				done
# 			fi
# 		fi
# 	done
# done


#### submit jobs ####

cd $OUTDIR

#submit some jobs at one time

# for pops in BERG BPNP
# for pops in CCRO CEBO 
# for pops in GRLL JENN KBPP
# for pops in KLDR MOSQ PRDF ROME SPVY SSSP

# for pops in KLDR MOSQ PRDF ROME SPVY SSSP 
# do
# 	cd $OUTDIR/jobs/
# 	for file in ${pops}*.fst.sh
# 	do
# 		cd $OUTDIR/errs
# 		sbatch $OUTDIR/jobs/$file
# 	done
# done

# Get results

count1=0
count2=0

mkdir fst/final_results
cd fst/final_results

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
