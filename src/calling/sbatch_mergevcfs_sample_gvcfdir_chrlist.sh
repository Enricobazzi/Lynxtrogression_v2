#!/bin/bash
#SBATCH --time=00-00:30:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20

module load picard

# input the sample name
sample=$1
# input gvcf directory
gvcfdir=$2
# input chromosome list
chromosomes=($(cat $3))

for chr in ${chromosomes[@]}; do
    # read in the gvcf files
    gvcfs=($(ls ${gvcfdir}/${sample}.${chr}.*.g.vcf))
    
    # run picard MergeVcfs
    java -jar $EBROOTPICARD/picard.jar MergeVcfs \
       $(for gvcf in ${gvcfs[@]}; do echo "I=${gvcf}"; done) \
       O=${gvcfdir}/${sample}.${chr}.g.vcf &
done

wait
