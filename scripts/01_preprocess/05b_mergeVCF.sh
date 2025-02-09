#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 150:00:00
#SBATCH --job-name=mergeVCF
#SBATCH -e mergeVCF
#SBATCH -o mergeVCF

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 07/10/22                  Last Modified: 08/24/22 ###
###########################################################################
###########################################################################
###                     mergeVCF.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR
#module load picard/2.18.17
module load htslib
module load samtools
module load bcftools

# Use BCFtools merge

#Options:
#        --force-samples                resolve duplicate sample names
#        --print-header                 print only the merged header and exit
#        --use-header <file>            use the provided header
#    -0  --missing-to-ref               assume genotypes at missing sites are 0/0
#    -f, --apply-filters <list>         require at least one of the listed FILTER strings (e.g. "PASS,.")
#    -F, --filter-logic <x|+>           remove filters if some input is PASS ("x"), or apply all filters ("+") [+]
#    -g, --gvcf <-|ref.fa>              merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,I16:sum,IDV:max,IMF:max
#    -i, --info-rules <tag:method,..>   rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [DP:sum,DP4:sum]
#    -l, --file-list <file>             read file names from the file
#    -m, --merge <string>               allow multiallelic records for <snps|indels|both|all|none|id>, see man page for details [both]
#        --no-version                   do not append version and command line to the header
#    -o, --output <file>                write output to a file [standard output]
#    -O, --output-type <b|u|z|v>        'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
#    -r, --regions <region>             restrict to comma-separated list of regions
#    -R, --regions-file <file>          restrict to regions listed in a file
#        --threads <int>                number of extra output compression threads [0]


### For PBrun ####
# Step0: zip individual vcfs and get index

#while read -a line
#do
#	bgzip -c ${line[0]}.raw_variants.vcf > ${line[0]}.raw_variants.vcf.gz
#	tabix -p vcf ${line[0]}.raw_variants.vcf.gz
#done < ../gatk/all202.sampleMap

# Doing it by chromosome (put all unplaced scaffolds in 1 file)

# for named chromosomes

while read -a chr
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --mem=170GB
#SBATCH -t 150:00:00
#SBATCH --job-name=${chr}_mergeVCF
#SBATCH -e ${chr}_mergeVCF
#SBATCH -o ${chr}_mergeVCF

cd $SLURM_SUBMIT_DIR
#module load picard/2.18.17
module load htslib
module load samtools
module load bcftools

cd /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/pbrun/gzipVCFs/

bcftools merge \
-l all202.sampleVCF.list \
-O z \
-r ${chr} \
--threads 40 \
-o ../chrVCFs/all202.raw_variants.${chr}.vcf" \
> /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_chr/${chr}_mergeVCF.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt


# Submit jobs #
while read -a chr
do
	cd /fs/ess/scratch/PAS1533/smathur/ibd/errors/per_chr/
	sbatch  /fs/ess/scratch/PAS1533/smathur/ibd/jobcodes/per_chr/${chr}_mergeVCF.sh
done < /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_named.txt

# For unplaced scaffolds
cd /fs/ess/scratch/PAS1533/smathur/ibd/variantCall/pbrun/gzipVCFs/

bcftools merge \
-l all202.sampleVCF.list \
-O z \
-R /fs/ess/scratch/PAS1533/smathur/ibd/lists/chr_unplaced.txt \
--threads 40 \
-o ../chrVCFs/all202.raw_variants.unplaced.vcf













