#!/bin/sh

module load cesga/2020
module load whatshap/1.1

# pop
pop=${1}
# chr
chr=${2}
# vcfdir
vcfdir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs
# i_vcf
i_vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.${chr}.vcf
# o_vcf
o_vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.${chr}.ps.vcf
# ref
ref=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_rufus_mLynRuf2.2/mLynRuf2.2.revcomp.scaffolds.fa
# bamlist
if [ $pop == 'lpa' ]; then
    bamlist=($(grep "lp_sm" data/bamlists/lp_ll_introgression.bamlist))
elif [ $pop == 'wel' ]; then
    bamlist=($(grep -E "ll_ki|ll_ur" data/bamlists/lp_ll_introgression.bamlist))
elif [ $pop == 'eel' ]; then
    bamlist=($(grep -E "ll_ya|ll_vl" data/bamlists/lp_ll_introgression.bamlist))
elif [ $pop == 'sel' ]; then
    bamlist=($(grep "ll_ca" data/bamlists/lp_ll_introgression.bamlist))
fi

whatshap phase \
    -o ${o_vcf} \
    --tag=PS \
    --reference=${ref} \
    ${i_vcf} \
    $(echo ${bamlist[@]})
