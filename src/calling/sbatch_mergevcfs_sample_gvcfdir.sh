#!/bin/bash
#SBATCH --time=00-01:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

module load picard

# input the sample name
sample=$1
# input gvcf directory
gvcfdir=$2

# read in the gvcf files
gvcfs=($(ls ${gvcfdir}/*.vcf | grep $sample))

# run picard MergeVcfs
java -jar $EBROOTPICARD/picard.jar MergeVcfs \
   $(for gvcf in ${gvcfs[@]}; do echo "I=${gvcf}"; done) \
   O=${gvcfdir}/${sample}.wholegenome.g.vcf.gz
