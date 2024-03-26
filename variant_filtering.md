## Variant Filtering

The variants contained in the VCF file were filtered following these criteria:

  1. Variants from repetitive and low complexity regions
  2. Indels and non-biallelic variants
  3. Substitutions from reference species (non variant SNPs with AF=1)
  4. Variant quality filters, as GATK standard practices
  5. Depth $$
  6. Missing Data $$

----

### Find repetitive and low complexity regions

To identify repeats and low complexity regions of the genome we used both [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) and [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/).

Repeats were identified and annotated by the people at CNAG when creating the reference genome using RepeatModeler. These are saved in the file `Repeats.4jb.gff3.gz`.

Low complexity regions of the genome were identified using RepeatMasker using the [mask_repeats](src/variant_filtering/mask_repeats.sh) script.
```
ref=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2/mLynRuf2.2.revcomp.scaffolds.fa

sbatch \
    --job-name=repeatmasker \
    --output=logs/variant_filtering/repeatmasker.out \
    --error=logs/variant_filtering/repeatmasker.err \
    src/variant_filtering/mask_repeats.sh ${ref}
```

This command generates the file `mLynRuf2.2.revcomp.scaffolds.fa.out.gff` which contains coordinates for low complexity regions.

To join the repeats with the low complexity regions:
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2

cat <(grep -v "#" ${ref_dir}/Repeats.4jb.gff3 | awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}') \
    <(grep -v "#" ${ref_dir}/repeatmasker/mLynRuf2.2.revcomp.scaffolds.fa.out.gff | awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}') |
    sort -k 1,1 -k2,2n -k3,3n |
    bedtools merge -i - \
    > ${ref_dir}/repeats_lowcomplexity_regions.bed
```

To calculate the length of these regions I run:
`awk '{sum+=$3-$2} END {print sum}' ${ref_dir}/repeats_lowcomplexity_regions.bed`

The result is 1_037_211_281, which is ~43% of the total genome (2_420_127_838).

----

### Applying filters 1 to 4

To apply the filters 1 to 4 I use the [apply_filters_1to4_invcf_ref_maskbed](src/variant_filtering/apply_filters_1to4_invcf_ref_maskbed.sh) script. This script applies a combination of [bedtools](https://bedtools.readthedocs.io/en/latest/), [gatk](https://gatk.broadinstitute.org/hc/en-us) and [bcftools](https://samtools.github.io/bcftools/bcftools.html) to remove:

  1. Variants from repetitive and low complexity regions
  2. Indels and non-biallelic variants
  3. Substitutions from reference species (non variant SNPs with AF=1)
  4. Variant quality filters, as GATK standard practices

```
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs
invcf=${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.vcf.gz
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
ref=${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa
maskbed=${ref_dir}/repeats_lowcomplexity_regions.bed

sbatch \
    --job-name=filter1to4 \
    --output=logs/variant_filtering/filter1to4.out \
    --error=logs/variant_filtering/filter1to4.err \
    src/variant_filtering/apply_filters_1to4_invcf_ref_maskbed.sh \
    ${invcf} ${ref} ${maskbed}
```
----

### Divide the VCF by population

The next filters are applied based on population so I will divide the filter4 vcf file into four population vcfs: `lpa`, `wel`, `eel` and `sel`.

I use gatk to do it:

```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
ref=${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs
invcf=${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.filter4.vcf

for pop in lpa wel eel sel; do
    if [ ${pop} == "lpa" ]; then
        samples=($(grep "lp_sm" data/sample.list))
    elif [ ${pop} == "wel" ]; then
        samples=($(grep -E "ll_ki|ll_ur" data/sample.list))
    elif [ ${pop} == "eel" ]; then
        samples=($(grep -E "ll_ya|ll_vl" data/sample.list))
    elif [ ${pop} == "sel" ]; then
        samples=($(grep -E "ll_ca" data/sample.list))
    fi
    echo "-- creating vcf of ${pop} --"
    gatk SelectVariants \
        -R ${ref} \
        -V ${invcf} \
        $(for sample in ${samples[@]}; do echo "-sn ${sample}";done) \
        -O ${invcf/.vcf/.${pop}_pop.vcf}
done
```

----

### Calculate missing data filters per population

To filter out excessively missing variants in each population I calculate missing data separately for each population.

----

### Calculate read depth filters in 10k bp window

To avoid including in the analysis possible paralogs whose SNP profiles do not reflect real genetic diversity we eliminate genomic regions where an excess of sequencing reads align to the reference genome.

Mean read depth in consecutive 10kbp windows along the genome was calculated using the software [mosdepth v0.3.2](https://github.com/brentp/mosdepth) from each sample's BAM using the [sbatch_mosdepth_10k_bam_outdir](src/variant_filtering/sbatch_mosdepth_10k_bam_outdir.sh) script:
```
inbams=($(cat data/bamlists/lp_ll_introgression.bamlist))

for bam in ${inbams[*]}; do
    echo "calculating depth of $bam"
    sbatch src/variant_filtering/sbatch_mosdepth_10k_bam_outdir.sh \
        ${bam} \
        data/variant_filtering/depth
done
```

From these results, we can calculate if a particular window should be filtered or not for a population using the [make_rdfilter_beds](src/variant_filtering/make_rdfilter_beds.py) scripts. Windows whose the sum of depth values for the population exceeds 1.5 times the mode of values are marked as failed for that population in the [regions_depth_filtering.tsv](data/variant_filtering/depth/regions_depth_filtering.tsv) table. The script also outputs one bed file for each population containing the coordinates for the windows that do not pass the filter, prints the following summary:
```
sel fail: 1050
wel fail: 873
eel fail: 1006
lpa fail: 1202
all fail: 714
lpa fail and others pass: 275
```
and these plots:

<head>
    <style>
        td, img {
            padding: 0;
            margin: 0;
        }
    </style>
</head>
<body>
    <table style="border-collapse: collapse; border-spacing: 0;">
        <tr>
            <td><img src="data/variant_filtering/depth/lpa_depth_distribution.png" alt="lpa_depth_distribution" style="width: 75%;" /></td>
            <td><img src="data/variant_filtering/depth/wel_depth_distribution.png" alt="wel_depth_distribution" style="width: 75%;" /></td>
        </tr>
        <tr>
            <td><img src="data/variant_filtering/depth/eel_depth_distribution.png" alt="eel_depth_distribution" style="width: 75%;" /></td>
            <td><img src="data/variant_filtering/depth/sel_depth_distribution.png" alt="sel_depth_distribution" style="width: 75%;" /></td>
        </tr>
    </table>
</body>
