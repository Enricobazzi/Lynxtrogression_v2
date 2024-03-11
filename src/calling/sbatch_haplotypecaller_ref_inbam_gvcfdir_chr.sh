#!/bin/bash
#SBATCH --time=00-06:00:00
#SBATCH --mem=60G
#SBATCH --output=logs/calling/gvcf-%j.out
#SBATCH --error=logs/calling/gvcf-%j.err
#SBATCH --cpus-per-task=24

# Purpose: call gvcf files for a bam file using HaplotypeCaller
# Prepare the bed files for each chromosome like this:
#  - cd /path/to/ref/folder
#  - mkdir chromosome_beds
#  - make a list of chromosome names in a file (e.g. chr.list)
#  - for chr in $(cat chr.list); do
#        grep ${chr} reference.fa.fai |
#            awk '{print $1, 0, $2}' | tr ' ' '\t' |
#            bedtools makewindows -b - -n 8 \
#            > chromosome_beds/${chr}.mk8win.bed
#        for n in {1..8}; do
#            sed "${n}q;d" chromosome_beds/${chr}.mk8win.bed \
#                > chromosome_beds/${chr}.${n}.bed
#        done
#    done
# The above code will create 8 bed files for each chromosome, each of which
# will contain 1/8th of the chromosome. This is to parallelize the gvcf
# calling process. The bed files will be used to call gvcf files for each
# 1/8th of the chromosome in parallel.

# Input:
#   - $1: reference genome
#   - $2: input bam file
#   - $3: output gvcf folder
#   - $4: chromosome name
# Output:
#   - gvcf files for each 1/8th of the chromosome
# Usage:
#   - bash call_gvcf_ref_bam_outgvcf_bed.sh path/to/ref.fa path/to/in.bam path/to/gvcf chr1

module load gatk

ref=$1
ref_folder=$(dirname ${ref})

inbam=$2
sample=$(basename -a ${inbam} | cut -d'_' -f1,2,3,4)

gvcf_dir=$3
chromosome=$4


for k in {1..8}; do
    out=${gvcf_dir}/${sample}.${chromosome}.${k}.g.vcf
    bed=${ref_folder}/chromosome_beds/${chromosome}.${k}.bed

    gatk HaplotypeCaller \
        -R ${ref} \
        -I ${inbam} \
        -L ${bed} \
        -O ${out} \
        --native-pair-hmm-threads 3 &
done

wait