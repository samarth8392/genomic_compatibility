#!/bin/sh -l
#SBATCH -A fnrtowhee
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=02_provean
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     02_provean.sh  				                    ###
###########################################################################

# see: https://github.com/cDaIryPOUV/ePat

# To get PROVEAN scores for the non-synonymous mutations in a dataset, we need three things.
# Reference genome (.fasta), genome annotation (.gff) and the SNPeff database you built for SNPeff (predictor.bin)
# A VCF file containing only non-synonymous mutations (i.e. missense mutations). 
#  ePAT installed in your working directory.


cd $SLURM_SUBMIT_DIR
# module load ondemand/project
# module load singularity/centos7
module load biocontainers
module load vcftools
module load snpeff
module load python

# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes"
SNPEFFDIR="/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/snpEff"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/provean"



# Step1: Divide all missense mutations into smaller sets

cd $OUTDIR
mkdir sites
cd sites
split -l 1000 -d $MAINDIR/sites/final152.missense.pos.txt snpSet

# Step2: Create folders for each file (N=69) with data and scripts
# NOTE: within each .vcf_dir, there should be the SNPeff data and bin and snpEff.jar file
# .
# ├── output
# │   └── output_provean_final152.nonSyn.snpSet00.xlsx
# └── snpEff
#     ├── build.logfile.txt
#     ├── data
#     │   └── scat
#     │       ├── genes.gff
#     │       ├── sequence.fa
#     │       └── snpEffectPredictor.bin
#     ├── final152.allchr.nonSyn.SNPEff.vcf
#     ├── final152.allchr.Syn.SNPEff.vcf
#     ├── no-downstream
#     ├── no-downstream.genes.txt
#     ├── provean_tmp
#     ├── snpEff.config
#     └── snpEff.jar


cd $OUTDIR
mkdir results

for i in {00..68}
do
	cd results/
	mkdir $i
	cd $i
	vcftools --vcf $SNPEFFDIR/final152.allchr.nonSyn.SNPEff.vcf \
	--recode --recode-INFO-all \
	--positions $OUTDIR/sites/snpSet$i \
	--out final152.nonSyn.snpSet$i

	mv final152.nonSyn.snpSet$i.recode.vcf final152.nonSyn.snpSet$i.vcf 
	cp -lr $OUTDIR/ePat/ ./

	mkdir final152.nonSyn.snpSet$i.vcf_dir
	cd final152.nonSyn.snpSet$i.vcf_dir

	cp -r $SNPEFFDIR ./
	cp $SNPEFFDIR/snpEff.config snpEff/

	cd $OUTDIR

 done

mkdir $OUTDIR/jobs
mkdir $OUTDIR/errors

for r in {00..68}
do
	echo "#!/bin/sh -l
#SBATCH -A fnrtowhee
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 150:00:00
#SBATCH --job-name=snp${r}.ePat
#SBATCH -e %x
#SBATCH -o %x


cd $SLURM_SUBMIT_DIR
# module load ondemand/project
# module load singularity/centos7

cd $OUTDIR/results/${r}/

i=final152.nonSyn.snpSet$i.vcf
SHARED_DIR=$OUTDIR/results/${r}/final152.nonSyn.snpSet${r}.vcf_dir/snpEff
WORK_DIR=$OUTDIR/results/${r}/
TMP_DIR=/tmp/final152.nonSyn.snpSet${r}.vcf
###############################################
mkdir -p /tmp/final152.nonSyn.snpSet${r}.vcf 

singularity run -B /tmp/final152.nonSyn.snpSet${r}.vcf:/root/tmp \
-B $OUTDIR/results/${r}/final152.nonSyn.snpSet${r}.vcf_dir/snpEff:/root/snpEff \
-B $OUTDIR/results/${r}/:$OUTDIR/results/${r}/ \
-W $OUTDIR/results/${r}/ \
$OUTDIR/results/${r}/ePat/ePat.sif \
$OUTDIR/results/${r}/ePat/script/automated_provean.sh \
-i final152.nonSyn.snpSet${r}.vcf -r scat 

rm -rf /tmp/final152.nonSyn.snpSet${r}.vcf" \
> $OUTDIR/jobs/ePat.snpSet${r}.sh

done

#### submit jobs ####

for i in {00..68}
do
	cd $OUTDIR/errors
	sbatch $OUTDIR/jobs/ePat.snpSet${i}.sh
done

# get results

cd $OUTDIR

for i in {00..68}
do
	cp results/$i/final152.nonSyn.snpSet${i}.vcf_dir/output/*xlsx final/xl/
	cp results/$i/final152.nonSyn.snpSet${i}.vcf_dir/output/*txt final/txt/
done

# Get final concatenated output

python $MAINDIR/jobcodes/scripts/python/concat_dataframe.py \
/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/provean/final/txt \
final152.nonSyn.provean

# Get deleterious and neutral functional variation

cut -f1,2,182,183 final152.nonSyn.provean.txt | grep "D" | sort -k1,1 -k2,2n > final152.deleterious.pos.txt

# Get the annotations for each nonSyn SNP

OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/provean"

cd $OUTDIR
cut -f1,2,8 final152.nonSyn.provean.txt > final152.nonSyn.annotate.txt


