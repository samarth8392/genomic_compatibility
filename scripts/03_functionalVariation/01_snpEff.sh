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
###########################################################################
###########################################################################
###                     01_snpEff.sh    		                    	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load java/21.0.2
module load vcftools

# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/snpeff"
SNPEFF="/fs/ess/scratch/PAS1533/smathur/genRes/annoVar/snpeff/snpEff"

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
# cp $SNPEFF/snpEff.config ./

# $SNPEFF/exec/snpeff build -noCheckCds -noCheckProtein \
# -c snpEff.config -gff3 -v scat &> build.logfile.txt

# If the database builds without error, you should see snpEffectPredictor.bin inside your scat folder

# Step6: Annotate your variants

$SNPEFF/exec/snpeff ann \
-stats -c snpEff.config \
-no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v scat \
$MAINDIR/vcf/final152.allchr.finalSNPs.vcf.gz \
> final152.allchr.SNPEff.vcf


# Step7: Use SNPSift to extract different synonymous and nonsynonymous SNPs
# see: https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files


cat final152.allchr.SNPEff.vcf | \
java -jar $SNPEFF/SnpSift.jar filter "( EFF[*].EFFECT = 'missense_variant' )" | \
grep -v "WARNING" > final152.allchr.missense.SNPEff.vcf

cat final152.allchr.SNPEff.vcf | \
java -jar $SNPEFF/SnpSift.jar filter "( EFF[*].EFFECT = 'synonymous_variant' )" | \
grep -v "WARNING" > final152.allchr.Syn.SNPEff.vcf

# high impact variants (LoF)
cat final152.allchr.SNPEff.vcf | \
java -jar $SNPEFF/SnpSift.jar filter "( EFF[*].IMPACT = 'HIGH' )" | \
grep -v "WARNING" > final152.allchr.nonsense.SNPEff.vcf


for type in nonsense Syn missense
do
	grep -v "#" final152.allchr.$type.SNPEff.vcf \
	| cut -f8 |  cut -f16 -d ";" | grep -v "WARNING" > $MAINDIR/sites/final152.$type.ann.txt

	grep -v "#" final152.allchr.$type.SNPEff.vcf \
	| grep -v "WARNING" | cut -f1,2  > $MAINDIR/sites/final152.$type.pos.txt

done

# Combine missense and nonsense into non-synonymous mutations

cd  $MAINDIR/sites/
cat final152.nonsense.pos.txt final152.missense.pos.txt > final152.nonSyn.pos.txt

# Get nonSyn vcf
cd $OUTDIR
vcftools --vcf final152.allchr.SNPEff.vcf \
--positions $MAINDIR/sites/final152.nonSyn.pos.txt \
--recode --recode-INFO-all \
--out final152.allchr.nonSyn.SNPEff