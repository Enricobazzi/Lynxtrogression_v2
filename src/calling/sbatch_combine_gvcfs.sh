#!/bin/bash
#SBATCH --time=03-00:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1

module load gatk

ref=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2/mLynRuf2.2.revcomp.scaffolds.fa
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_gvcfs
gvcfs=($(ls ${gvcf_dir}/*.wholegenome.g.vcf.gz))
out_gvcf=${gvcf_dir}/lynxtrogression_v2.g.vcf.gz

gatk CombineGVCFs \
    -R ${ref} \
    $(for gvcf in ${gvcfs[@]}; do echo "--variant ${gvcf}"; done) \
    -O ${out_gvcf}
