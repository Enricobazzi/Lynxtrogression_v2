#!/bin/bash
#SBATCH --time=01-00:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1

module load gatk

# input the gvcf list
gvcfs=($(cat $1))
# input the reference genome
ref=$2
# input the output gvcf
out_gvcf=$3

gatk CombineGVCFs \
    -R ${ref} \
    $(for gvcf in ${gvcfs[@]}; do echo "--variant ${gvcf}"; done) \
    -O ${out_gvcf}
