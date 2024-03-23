#!/bin/bash
#SBATCH --time=01-00:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

module load gatk

# read the reference genome
ref=$1
# read the input gvcf
ingvcf=$2
# read the output vcf
outvcf=$3

# run GATK GenotypeGVCFs
gatk GenotypeGVCFs \
    -R ${ref} \
    -V ${ingvcf} \
    -O ${outvcf}
