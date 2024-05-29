#!/bin/bash
#SBATCH -c 10
#SBATCH --mem=40G
#SBATCH -t 03-00:00:00

# load the necessary modules
module load cesga/system miniconda3/22.11.1-1
source activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/gadma2

cd data/demographic_inference

# define the pair, the iteration number and the folder
pair=$1
n=$2
folder=$3

# run gadma
if [ -d ${folder} ]; then
    gadma --resume ${folder}
else
    gadma -p ${pair}.gadma.parameters -o ${pair}_${n}
fi
