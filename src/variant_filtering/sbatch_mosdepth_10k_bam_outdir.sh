#!/bin/bash
#SBATCH -e logs/variant_filtering/mosdepth-%j.err
#SBATCH -o logs/variant_filtering/mosdepth-%j.out
#SBATCH --time=00-00:30:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10

module load mosdepth

# input bam file
bam=$1

# sample name
sample=$(basename -a ${bam} | cut -d'_' -f1,2,3,4)

# output directory
outdir=$2

# run mosdepth
mosdepth -t 10 -n --by 10000 ${outdir}/${sample} ${bam}
