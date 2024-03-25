#!/bin/bash
#SBATCH --time=00-06:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

module load bedtools
module load gatk
module load bcftools

# read the input vcf
invcf=$1
# read the reference genome
ref=$2
# read bed file with repetitive and low complexity regions
maskbed=$3

# declare the output vcf from the input vcf
outvcf=$(dirname ${invcf})/$(basename ${invcf%.vcf.gz})

# apply filter 1: remove variants from repetitive and low complexity regions
bedtools subtract -a ${invcf} -b ${maskbed} -header | uniq > ${outvcf}.filter1.vcf

# apply filter 2: remove INDELs and non-biallelic variants
gatk SelectVariants \
    -R ${ref} \
    -V ${outvcf}.filter1.vcf \
    -O ${outvcf}.filter2.vcf \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC

# apply filter 3: remove variants with AF=1
bcftools view -e 'INFO/AF=1.00' ${outvcf}.filter2.vcf > ${outvcf}.filter3.vcf

# apply filter 4: remove low quality variants (as defined by GATK)
bcftools view -e 'INFO/MQRankSum < -12.5 | INFO/ReadPosRankSum < -8.0 | QUAL < 30 | INFO/QD < 2.0 | INFO/FS > 60.0 | INFO/MQ < 40.0' \
    ${outvcf}.filter3.vcf > ${outvcf}.filter4.vcf
