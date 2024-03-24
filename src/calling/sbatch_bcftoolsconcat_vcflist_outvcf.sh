#!/bin/bash
#SBATCH --time=00-00:10:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=22

# load bcftools
module load bcftools

# read the list of vcfs
vcflist=$1
# read the output vcf
outvcf=$2

# run bcftools concat
bcftools concat \
    -f ${vcflist} \
    -O z -o ${outvcf} \
    --threads 22