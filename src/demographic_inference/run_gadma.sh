#!/bin/bash
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH -t 3-00:00:00

# load the necessary modules
module load cesga/system miniconda3/22.11.1-1
source activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/gadma2

cd data/demographic_inference

# define the pair and the iteration number
pair=$1
n=$2

# run gadma
gadma -p ${pair}.gadma.parameters -o ${pair}_${n}
