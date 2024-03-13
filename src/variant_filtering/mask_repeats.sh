#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=12

module load repeatmasker

# read reference genome as first argument
ref=$1

# dir is created in the same directory as the reference genome
ref_dir=$(dirname $ref)
if [ ! -d ${ref_dir}/repeatmasker ]; then
    mkdir ${ref_dir}/repeatmasker
fi

RepeatMasker -pa 12 \
    -gff -xsmall -noint \
    -dir ${ref_dir}/repeatmasker \
    ${ref}
