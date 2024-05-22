## Demographic Inference

I will be using [GADMA2](https://github.com/ctlab/GADMA) (see [Noskova et al. 2023](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giad059/7248629)) to reconstruct the demographic history of the Iberian and Eurasian lynx population pairs.

----

### Prepare the Dataset - Remove genes

I will use neutral regions of the genome and independent SNPs for more solid reconstruction of the neutral demographic history

I extract the genes as a BED file from the GFF annotation of the [Lynx rufus reference genome](https://denovo.cnag.cat/lynx_rufus)

```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2

awk '$3 == "gene" {print}' ${ref_dir}/LYRU2_2A.FA.gff3 |
    cut -f1,4,5,9 | awk '{print $1, $2-1, $3, $4}' | tr ' ' '\t' \
    > ${ref_dir}/genes.bed
```

Then I subtract them from the VCF:
```
genes_bed=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2/genes.bed
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

for pair in lpa-wel lpa-eel lpa-sel; do
    echo "removing genes from ${pair}_pair.miss_fil.rd_fil.variant vcf"
    bedtools subtract -header \
        -a ${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pair}_pair.miss_fil.rd_fil.variant.vcf \
        -b ${genes_bed} \
    > ${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pair}_pair.miss_fil.rd_fil.variant.nogenes.vcf
done
```

----

### Prepare the Dataset - Prune SNPs

I prune the SNPs using `vcftools`, keeping SNPs with no missing data and at a distance of at least 50kb between eachother:
```
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

for pair in lpa-wel lpa-eel lpa-sel; do
    vcf=${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pair}_pair.miss_fil.rd_fil.variant.nogenes.vcf
    
    echo "pruning ${pair}"
    vcftools --vcf ${vcf} \
        --thin 50000 --max-missing 1 \
        --recode --recode-INFO-all \
        --out data/demographic_inference/${pair}.pruned
done
```

This leaves the following amount of SNPs in each population pair (info in data/demographic_inference/ logs):
```
# lpa-wel: 31628 out of 2735580
# lpa-eel: 31885 out of 3004621
# lpa-sel: 31882 out of 2732135
```
----

### Prepare the Dataset - Calculate Observed Genome Length

To calculate the amount of observed genome I need to subtract from the autosome length all of the regions I removed during filtering (low complexity and repetitive, high depth, genes) and also keep a proportion of sites equal to the proportion of SNPs I kept when filtering out SNPs (low quality, missing, pruned).

To calculate these:
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2

for pop in wel eel sel; do
    echo "lpa-${pop}:"
    bedtools subtract \
        -a <(grep -i "chr" ${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa.fai | grep -viE "chry|chrx" | awk '{print $1, 0, $2}' | tr ' ' '\t') \
        -b ${ref_dir}/repeats_lowcomplexity_regions.bed |
    bedtools subtract \
        -a stdin \
        -b ${ref_dir}/genes.bed |
    bedtools subtract \
        -a stdin \
        -b <(bedtools merge -i <(cat data/variant_filtering/depth/lpa.rd_filter.bed data/variant_filtering/depth/${pop}.rd_filter.bed | sort -k1,1 -k2,2n)) \
    > data/demographic_inference/lpa-${pop}.callable.bed
    cat data/demographic_inference/lpa-${pop}.callable.bed | awk '{sum += $3 - $2} END {print sum}'
done
```

So the total genomic region observed is:
```
# lpa-wel: 735474363
# lpa-eel: 735395197
# lpa-sel: 735147869
```

From these regions I take away a proportion of SNPs for low quality, a proportion for missing data and a proportion for pruning:
```
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

# proportion of low quality variants: 6495225 / 7131855 = 0.910734304048526
grep -v "#" ${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter3.vcf | wc -l
# 7131855
grep -v "#" ${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.vcf | wc -l
# 6495225

# proportion of highly missing variants in each pop-pair
# lpa-wel: 6117961 / 6495225 = 0.9419167157411791
# lpa-eel: 6113904 / 6495225 = 0.9412921030449292
# lpa-sel: 6198505 / 6495225 = 0.9543172099503866

for pair in lpa-wel lpa-eel lpa-sel; do
    echo "lpa-${pop}:"
    grep -v "#" ${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pair}_pair.vcf | wc -l
    grep -v "#" ${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pair}_pair.miss_fil.vcf | wc -l
done

# proportion of pruned SNPs
# lpa-wel: 31628 / 2735580 = 0.011561716345345412
# lpa-eel: 31885 / 3004621 = 0.01061198733550754
# lpa-sel: 31882 / 2732135 = 0.011669262316832806

for pair in lpa-wel lpa-eel lpa-sel; do
    echo "lpa-${pop}:"
    grep -v "#" ${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pair}_pair.miss_fil.rd_fil.variant.nogenes.vcf | wc -l
    grep -v "#" data/demographic_inference/${pair}.pruned.recode.vcf | wc -l
done
```

The total number of observed sites in each population pair is:
```
# lpa-wel: 735474363 * 0.910734304048526 * 0.9419167157411791 * 0.011561716345345412 = 7294475.137109491
# lpa-eel: 735395197 * 0.910734304048526 * 0.9412921030449292 * 0.01061198733550754 = 6690115.6057525985
# lpa-sel: 735147869 * 0.910734304048526 * 0.9543172099503866 * 0.011669262316832806 = 7455942.606560742
```

----

### Installing GADMA2

I create a conda environment for easier GADMA installation. As dependecies are tricky with the installation, I found the following to work:
```
cd /mnt/netapp1/Store_CSIC/home/csic/eye/eba # $STORE

module load cesga/system miniconda3/22.11.1-1

conda create --prefix=gadma2 python=3.10
source activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/gadma2

pip install git+https://github.com/MomentsLD/moments.git
conda install -c conda-forge dadi
conda install -c conda-forge scikit-allel
pip install -i https://test.pypi.org/simple/ gadma
pip uninstall ruamel.yaml
pip install "ruamel.yaml<0.18.0"
pip uninstall matplotlib
pip install "matplotlib<3.5"
```

where `pip install -i https://test.pypi.org/simple/ gadma` installs a development version of GADMA where `Lower bound of first split` can be defined. This is ideal with my data, two populations of two distinct species, where the time intervals for within-population events can be very small (e.g. recent population decline), but the divergence time might be very high.

I can activate this environment by running the following:
```
module load cesga/system miniconda3/22.11.1-1 && source activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/gadma2
```

### Preps

To run GADMA I need 3 things:
- vcf file
- popmap file
- parameter file

Run gadma:
```
pair=lpa-wel
for n in {1..50}; do
    sbatch --job-name=${n}_gadma_${pair} \
        --output=logs/demographic_inference/${pair}.gadma.${n}.out \
        --error=logs/demographic_inference/${pair}.gadma.${n}.err \
        src/demographic_inference/run_gadma.sh ${pair} ${n}
done
```