## Demographic Inference

I will be using [GADMA2](https://github.com/ctlab/GADMA) (see [Noskova et al. 2023](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giad059/7248629)) to reconstruct the demographic history of the Iberian and Eurasian lynx population pairs.


### Prepare the Dataset - Remove genes & Prune SNPs

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

### Prune SNPs

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
- lpa-wel: 31628 out of 2735580
- lpa-eel: 31885 out of 3004621
- lpa-sel: 31882 out of 2732135

### Calculate Observed Genome Length


### Run GADMA2

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
