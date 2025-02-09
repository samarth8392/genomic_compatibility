#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 150:00:00
#SBATCH --job-name=snpEff
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/02/22                  Last Modified: 09/22/23 ###
###########################################################################
###########################################################################
###                     01_snpEff.sh    		                    	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load snpeff
module load vcftools

# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/snpeff"
SNPEFF="/usr/local/snpeff/4.2/snpEff"

# SNPeff

# Steps:
# 1. Copy snpeff.config file from /group/bioinfo/apps/apps/snpEff-4.3/ to local directory
# 2. Make a directory in your local directory called "data" and another directory within data with your species name. I call it scat (S.catanatus)
# 3. copy ref fasta and annotation gff as sequence.fa and genes.gff into "scat" folder
# 4. Change the following in config file:
# 	(a) Add in the first two lines:
#			# S.catanatus genome
#			scat.genome : S.catanatus
#	(b) Change data.dir to $OUTDIR/data
# 5. Build the database

cd $OUTDIR
#cp $SNPEFF/snpEff.config ./

# java -jar $SNPEFF/snpEff.jar build \
# -c snpEff.config -gff3 -v scat &> build.logfile.txt

# If the database builds without error, you should see snpEffectPredictor.bin inside your scat folder

# Step6: Annotate your variants

# java -jar $SNPEFF/snpEff.jar ann \
# -stats -c snpEff.config \
# -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v scat \
# $MAINDIR/vcf/final172.allchr.finalSNPs.vcf.gz \
# > final172.allchr.SNPEff.vcf


# Step7: Use SNPSift to extract different synonymous and nonsynonymous SNPs
# see: https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files


# cat final172.allchr.SNPEff.vcf | \
# java -jar $SNPEFF/SnpSift.jar filter "( EFF[*].EFFECT = 'missense_variant' )" | \
# grep -v "WARNING" > final172.allchr.missense.SNPEff.vcf

# cat final172.allchr.SNPEff.vcf | \
# java -jar $SNPEFF/SnpSift.jar filter "( EFF[*].EFFECT = 'synonymous_variant' )" | \
# grep -v "WARNING" > final172.allchr.Syn.SNPEff.vcf

# # high impact variants (LoF)
# cat final172.allchr.SNPEff.vcf | \
# java -jar $SNPEFF/SnpSift.jar filter "( EFF[*].IMPACT = 'HIGH' )" | \
# grep -v "WARNING" > final172.allchr.nonsense.SNPEff.vcf


# for type in nonsense Syn missense
# do
# 	grep -v "#" final172.allchr.$type.SNPEff.vcf \
# 	| cut -f8 |  cut -f16 -d ";" | grep -v "WARNING" > $MAINDIR/sites/final172.$type.ann.txt

# 	grep -v "#" final172.allchr.$type.SNPEff.vcf \
# 	| grep -v "WARNING" | cut -f1,2  > $MAINDIR/sites/final172.$type.pos.txt

# done

# Combine missense and nonsense into non-synonymous mutations

# cd  $MAINDIR/sites/
# cat final172.nonsense.pos.txt final172.missense.pos.txt > final172.nonSyn.pos.txt

# Get nonSyn vcf
cd $OUTDIR
vcftools --vcf final172.allchr.SNPEff.vcf \
--positions $MAINDIR/sites/final172.nonSyn.pos.txt \
--recode --recode-INFO-all \
--out final172.allchr.nonSyn.SNPEff

# After filtering, kept 69,262 out of a possible 18,923,075 Sites

### SYNONYMOUS = 93,195
### MISSENSE = 68,041
### NONSENSE = 1,221
### NONSYNONYMOUS = 69,262