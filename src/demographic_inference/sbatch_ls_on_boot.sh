#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=60G
#SBATCH -t 00-06:00:00

# load the necessary modules
module load cesga/system miniconda3/22.11.1-1
source activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/gadma2

# read the arguments
boot_folder=$1
best_moments=$2
output_folder=$3

if [ ! -d ${output_folder} ]; then
    mkdir -p ${output_folder}
fi

gadma-run_ls_on_boot_data \
    -b ${boot_folder} \
    -d ${best_moments} \
    -o ${output_folder} \
    -j 20
