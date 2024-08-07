#!/bin/sh
module load cesga/2020
module load gcccore/system shapeit4/4.2.1

# pop
pop=${1}
# chr
chr=${2}
# vcfdir
vcfdir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

# i_vcf
i_vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.${chr}.ps.vcf.gz
# gmap
gmap=data/phasing/${chr}.gmap
# o_vcf
o_vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.${chr}.ps.phased.vcf

shapeit4.2 \
 --input ${i_vcf} \
 --map ${gmap} \
 --region ${chr} \
 --use-PS 0.0001 \
 --output ${o_vcf} \
 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
