#!/bin/bash
#SBATCH -o logs/fastp/fastp-%j.out
#SBATCH -e logs/fastp/fastp-%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=6

module load cesga/2020 gcccore/system fastp/0.22.0

# This script takes three positional arguments:
# 1. the directory where the fastq files are located
# 2. the r1 fastq file
# 3. the r2 fastq file
# and runs fastp to trim adapters and low quality reads

# read the fastq directory from the first positional argument
fastq_dir=${1}
# read the r1 from the second positional argument
r1_fastq=${2}
# read the r2 from the third positional argument
r2_fastq=${3}

# make fastq files have the same extension .fastp.fastq.gz 
# regardless of the original extension .fastq.gz or .fq.gz
# if the original is neither then the script will fail
if [[ ${r1_fastq} == *.fastq.gz ]]; then
    r1_fastp=${r1_fastq/.fastq.gz/.fastp.fastq.gz}
    r2_fastp=${r2_fastq/.fastq.gz/.fastp.fastq.gz}
elif [[ ${r1_fastq} == *.fq.gz ]]; then
    r1_fastp=${r1_fastq/.fq.gz/.fastp.fastq.gz}
    r2_fastp=${r2_fastq/.fq.gz/.fastp.fastq.gz}
else
    echo "Error: fastq files should have extension .fastq.gz or .fq.gz"
    exit 1
fi

# create the fastp directory if it does not exist
mkdir -p ${fastq_dir}/fastp

# run fastp
fastp \
    -i ${fastq_dir}/${r1_fastq} -I ${fastq_dir}/${r2_fastq} \
    -o ${fastq_dir}/fastp/${r1_fastp} -O ${fastq_dir}/fastp/${r2_fastp} \
    -h ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_fastp.html} -j ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_fastp.json} \
    --unpaired1 ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_unpaired.fastq.gz} --unpaired2 ${fastq_dir}/fastp/${r2_fastq/.fastq.gz/_unpaired.fastq.gz} \
    --failed_out ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_failed.fastq.gz} \
    --dont_overwrite \
    --trim_poly_g \
    --length_required 30 \
    --correction \
    --detect_adapter_for_pe \
    --thread 6
