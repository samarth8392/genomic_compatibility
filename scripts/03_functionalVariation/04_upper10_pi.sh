#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=upper10
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     04_upper10.sh              						###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load vcftools
module load htslib
module load samtools

# dirs
MAINDIR="/fs/ess/scratch/PAS1533/smathur/genRes"
OUTDIR="/fs/ess/scratch/PAS1533/smathur/genRes/upper10"
REFDIR="/fs/ess/scratch/PAS1533/smathur/genRes/ref/"


###################### 
# ADAPTIVE DIVERSITY #
###################### 

# Upper10% Pi Using ANGSD

# Step0a: Make BED of all upper10% genes

cd $OUTDIR
awk '{gsub(/;$/,"");print}' $REFDIR/annotations/Scate_HiC_rnd4.all.putative.function.gff > Scate_HiC.revisedgenome.all.gff
cat Scate_HiC.revisedgenome.all.gff | awk '$3 =="gene" {print $0}' > Scate_HiC.revisedgenome.all_justgenes.gff
cat Scate_HiC.revisedgenome.all.gff | awk '$3 =="CDS" {print $0}' > Scate_HiC.revisedgenome.all_justCDS.gff
cut -f1,4,5 Scate_HiC.revisedgenome.all_justgenes.gff > Scate_HiC.justGenes.bed
cut -f1,4,5 Scate_HiC.revisedgenome.all_justCDS.gff > Scate_HiC.justCDS.bed

# Step0b: Create bam lists

while read -a pop
do
    infile="$MAINDIR/lists/byPop/$pop.final152.sample.list"
    outfile1="$MAINDIR/lists/bamlists/$pop.final152.bamlist"
    outfile2="$MAINDIR/lists/bamlists/$pop.up10.bamlist"
    bamDir1="$MAINDIR/bams/final_bam"
    bamDir2="$MAINDIR/bams/upper10"

    touch $outfile1
    touch $outfile2

    while read -a sample
    do
        echo "$bamDir1/${sample}_mark_dups.bam" >> $outfile1
        echo "$bamDir2/${sample}.in.upper10.genes.bam" >> $outfile2
    done < $infile
done < $MAINDIR/lists/final152.onlyScat.popNames.txt

# Step0c: Reindex Fastas

cd $REFDIR
samtools faidx Scatenatus_HiC_v1.1.fasta


# Step1: Filter BAMs to retain only upper10% genes

cd /fs/ess/scratch/PAS1533/smathur/genRes/bams/final_bam

for sample in $(ls *bam| cut -f1 -d "_")
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 150:00:00
#SBATCH --job-name=${sample}.bamFilt
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge
module load samtools
module load htslib

cd /fs/ess/scratch/PAS1533/smathur/genRes/bams/upper10

samtools view /fs/ess/scratch/PAS1533/smathur/genRes/bams/final_bam/${sample}_mark_dups.bam \
-b -h -o ${sample}.in.upper10.genes.bam -U ${sample}.out.upper10.genes.bam \
-L $OUTDIR/beds/upp10_geneCoord.bed" \
>$OUTDIR/${sample}.up10.bamFilter.sh
done

Submit jobs
cd /fs/ess/scratch/PAS1533/smathur/genRes/bams/final_bam
for sample in $(ls *bam| cut -f1 -d "_")
do
	cd $OUTDIR/errs/
	sbatch  $OUTDIR/${sample}.up10.bamFilter.sh
done

index bam files
cd /fs/ess/scratch/PAS1533/smathur/genRes/bams/upper10/
for sample in $(ls *in*bam| cut -f1 -d ".")
do
	samtools index -b ${sample}.in.upper10.genes.bam
done


# index sites

cd $MAINDIR/sites/
~/softwares/angsd/angsd sites index final152.missense.pos.txt

# Step2: Calculate theta for each population in just the CDS of upper10% genes using ANGSD
# Step3: Calculate theta for each population in just the missense of upper10% genes using ANGSD

while read -a pop
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 150:00:00
#SBATCH --job-name=${pop}.up10
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge
module load samtools
module load gcc-compatibility
module load htslib

cd $OUTDIR/pi

# Get 1D SFS

mkdir ${pop}

~/softwares/angsd/angsd -P 128 \
-dosaf 1 -GL 1 \
-anc $REFDIR/Scatenatus_HiC_v1.1.fasta \
-bam $MAINDIR/lists/bamlists/$pop.up10.bamlist \
-rf $OUTDIR/beds/Scate_HiC.justCDS.regionsList.txt \
-out ${pop}/${pop}.upper10_justCDS 

~/softwares/angsd/misc/realSFS -P 128 -fold 1 \
${pop}/${pop}.upper10_justCDS.saf.idx -maxIter 100 \
> ${pop}/${pop}.upper10_justCDS.sfs

# Get thetas

~/softwares/angsd/misc/realSFS saf2theta \
${pop}/${pop}.upper10_justCDS.saf.idx \
-sfs ${pop}/${pop}.upper10_justCDS.sfs -fold 1 \
-outname ${pop}/${pop}.upper10_justCDS

~/softwares/angsd/misc/thetaStat do_stat \
${pop}/${pop}.upper10_justCDS.thetas.idx \
-outnames ${pop}/${pop}.upper10_justCDS

~/softwares/angsd/misc/thetaStat print \
${pop}/${pop}.upper10_justCDS.thetas.idx \
> ${pop}/${pop}.upper10_justCDS_persite.txt

# Dump Counts
~/softwares/angsd/angsd -P 128 \
-doCounts 1 -dumpCounts 2 \
-rf $OUTDIR/beds/Scate_HiC.justCDS.regionsList.txt \
-bam $MAINDIR/lists/bamlists/$pop.up10.bamlist  \
-out ${pop}/${pop}.upper10_justCDS


### Within noSyn sites

~/softwares/angsd/angsd -P 128 \
-dosaf 1 -GL 1 \
-anc $REFDIR/Scatenatus_HiC_v1.1.fasta \
-bam $MAINDIR/lists/bamlists/$pop.up10.bamlist \
-sites $MAINDIR/sites/final152.missense.pos.txt \
-out ${pop}/${pop}.upper10_missense

~/softwares/angsd/misc/realSFS -P 128 -fold 1 \
${pop}/${pop}.upper10_missense.saf.idx -maxIter 100 \
> ${pop}/${pop}.upper10_missense.sfs

# Get thetas

~/softwares/angsd/misc/realSFS saf2theta \
${pop}/${pop}.upper10_missense.saf.idx \
-sfs ${pop}/${pop}.upper10_missense.sfs -fold 1 \
-outname ${pop}/${pop}.upper10_missense

~/softwares/angsd/misc/thetaStat do_stat \
${pop}/${pop}.upper10_missense.thetas.idx \
-outnames ${pop}/${pop}.upper10_missense

~/softwares/angsd/misc/thetaStat print \
${pop}/${pop}.upper10_missense.thetas.idx \
> ${pop}/${pop}.upper10_missense_persite.txt

# Dump Counts
~/softwares/angsd/angsd -P 128 \
-doCounts 1 -dumpCounts 2 \
-sites $MAINDIR/sites/final152.missense.pos.txt \
-bam $MAINDIR/lists/bamlists/$pop.up10.bamlist \
-out ${pop}/${pop}.upper10_missense_sites" \
> $OUTDIR/jobs/${pop}.upper10.theta.sh
done < $MAINDIR/lists/final152.onlyScat.popNames.txt

# Submit jobs
while read -a pop
do
	cd $OUTDIR/errs/
	sbatch $OUTDIR/jobs/${pop}.upper10.theta.sh
done < $MAINDIR/lists/final152.onlyScat.popNames.txt



# Grab results

cd $OUTDIR/pi/summary

while read -a pop
do
	cp ../$pop/*.pestPG ./
	sed -i '1d' $pop.upper10_justCDS.pestPG
	sed -i '1d' $pop.upper10_missense.pestPG

done < $MAINDIR/lists/final152.onlyScat.popNames.txt


