# Phasing variants

Phasing of variants will be conducted with a pipeline that first uses [WhatsHap v.1.1](https://whatshap.readthedocs.io/en/latest/index.html) ([Martin et al., 2016](https://www.biorxiv.org/content/10.1101/085050v2)) to create phase sets from individual read and population data. The output of WhatsHap is then passed to [SHAPEIT4](https://odelaneau.github.io/shapeit4/) v.4.2.1 ([Delaneau et al., 2019](https://www.nature.com/articles/s41467-019-13225-y)) that will infer the haplotypes of each sample for each chromosome.

All this is based on what [Lorena already ran]((https://github.com/lorenalorenzo/PlanNacional_Selectionscans/blob/main/Data%20preparation/02_phasing.md)) for the samples mapped to the Felix catus reference genome.

### Splitting the VCF

To divide my VCF into single population VCFs and further dividing those into single chromosome VCFs I use bedtools:
```
# vcfdir
vcfdir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs
# refdir
refdir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))

# iterate pops
for pop in lpa wel eel sel; do
    vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.vcf
    for chr in ${chrs[@]}; do
        echo "extracting ${chr} from:"
        echo "$vcf"
        chrbed=${refdir}/chromosome_beds/${chr}.bed
        bedtools intersect -header -a ${vcf} \
            -b ${chrbed} \
        > ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.${chr}.vcf
    done
done
```

### Generate genetic map

To run SHAPEIT4 I also need to provide a genetic map for the SNPs to phase. As we don't have one, we will manually generate a genetic map by multiplying the physical distance in bp between SNPs and genome wide average recombination rate, which is 1.9 cM/Mbp. By cumulatively summing the multiplication of the physical distance from previous the SNP by 0.0000019, we obtain the cM value of each SNP. This approximation is not ideal but it's the only way we can provide a map. To calculate this I run:
```
# vcfdir
vcfdir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs
# refdir
refdir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))
# vcf
vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.vcf

for chr in ${chrs[@]}; do
    echo "calculating genetic map of ${chr} from:"
    echo "${vcf}"
    grep -v "#" ${vcf} | grep -w ${chr} | cut -f1,2 |
        awk '{ print $2, $1 }' |
        awk {'if ( NR==1 ) print $1, $2, 0; else print $1, $2, $1-p, ($1-p)*0.0000019; p=$1'} |
        awk 'BEGIN{print "pos", "chr", "cM"} {sum+=$4} {print $1, $2, sum}' |
        tr ' ' '\t' > data/phasing/${chr}.gmap
done
```
which will output a gmap table for each chromosome, made of 3 columns: position, chromosome, cM (format useful for SHAPEIT4).

### Generate phase sets with WhatsHap

For more precise phasing, we first run the software WhatsHap using the --tag=PS (see [here](https://whatshap.readthedocs.io/en/latest/guide.html#representation-of-phasing-information-in-vcfs)).

Phase sets are generated from the VCF of each chromosome of each population using the [run_whatshap_ps.sh](src/phasing/run_whatshap_ps.sh) script:
```
# refdir
refdir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))

# for pop in lpa wel sel eel; do
# for pop in lpa wel; do
for pop in sel eel; do
    for chr in ${chrs[@]}; do
        sbatch \
            --job-name=${pop}_${chr}_ps \
            --output=logs/phasing/${pop}_${chr}_ps.out \
            --error=logs/phasing/${pop}_${chr}_ps.err \
            -c 2 --mem 20G -t 03-00:00:00 \
            src/phasing/run_whatshap_ps.sh ${pop} ${chr}
    done
done
```

### Phase using SHAPEIT4

SHAPEIT needs the VCFs zipped and indexed:

```
# refdir
refdir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))
# vcfdir
vcfdir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

# for pop in lpa wel sel eel; do
# for pop in lpa wel; do
for pop in sel eel; do
    for chr in ${chrs[@]}; do
        # i_vcf
        i_vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.${chr}.ps.vcf
        echo "zipping ${i_vcf}"
        bgzip ${i_vcf}
        echo "indexing ${i_vcf}.gz"
        bcftools index ${i_vcf}.gz
    done
done
```
Run SHAPEIT4 to phase the gzipped, indexed, phase set vcf using the [run_shapeit.sh](src/phasing/run_shapeit.sh) script:
```
# refdir
refdir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))

# for pop in lpa wel sel eel; do
# for pop in sel eel; do
for pop in lpa wel; do
    for chr in ${chrs[@]}; do
        sbatch \
            --job-name=${pop}_${chr}_shapeit \
            --output=logs/phasing/${pop}_${chr}_shapeit.out \
            --error=logs/phasing/${pop}_${chr}_shapeit.err \
            -c 2 --mem 10G -t 00-01:00:00 \
            src/phasing/run_shapeit.sh ${pop} ${chr}
    done
done
```
Make pop_pair vcf by merging and concatenating chromosome vcf of each pop:
```
module load cesga/2020 bcftools/1.19

# refdir
refdir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))
# vcfdir
vcfdir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

## zip & index pops
# for pop in lpa wel sel eel; do
# for pop in sel eel; do
for pop in lpa wel; do
    for chr in ${chrs[@]}; do
        vcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.${pop}_pop.${chr}.ps.phased.vcf
        echo "zip & index: ${vcf}"
        bgzip ${vcf}
        tabix -p vcf ${vcf}.gz
    done
done

## merge pop pair
# for pop in wel eel sel
for pop in wel; do
    for chr in ${chrs[@]}; do
        echo "merging and indexing ${pop_pair} ${chr} vcf:"
        echo "${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.${chr}.ps.phased.merged.vcf"
        # merge chromosome vcf using bcftools (tac to reverse order so that lpa is pop2)
        bcftools merge -O "z" $(ls ${vcfdir}/*.ps.phased.vcf.gz | grep -E "${pop}_pop|lpa_pop" | grep ".${chr}." | tac) \
            > ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.${chr}.ps.phased.merged.vcf.gz
        # index vcf
        tabix -p vcf ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.${chr}.ps.phased.merged.vcf.gz
    done
done

## concat pop pair
# for pop in wel eel sel
for pop in wel; do
    echo "concatenating lpa-${pop}"
    vcf-concat $(ls ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.*.ps.phased.merged.vcf.gz) \
        > ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.vcf
done

## add old header to vcf
# for pop in wel eel sel
for pop in wel; do
    grep "##" ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.vcf \
        > ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.vcf
    echo "##INFO=<ID=CM,Number=1,Type=Float,Description='CentiMorgan'>" \
        >> ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.vcf
    grep -m1 "#CHR" ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.vcf \
        >> ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.vcf
    grep -v "#" ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.vcf \
        >> ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.vcf
done

## annotate AF,AN fields
# for pop in wel eel sel
for pop in wel; do
    echo "annotating lpa-${pop}"
    bcftools +fill-tags ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.vcf -- -t AF,AN \
    > ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.vcf
done
```

### Filter phased VCF

Apply same filters as other VCF (depth & non-variant sites - missing data removed by phasing):

```
module load cesga/2020 bcftools/1.19
module load bedtools

# vcfdir
vcfdir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

# for pop_pair in wel eel sel
for pop in wel; do
    # apply the filter based on read depth
    echo "filtering lpa-${pop} vcf removing windows with an excess of read depth:"
    bedtools subtract -header \
        -a ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.vcf \
        -b data/variant_filtering/depth/lpa.rd_filter.bed |
    bedtools subtract -header \
        -a stdin \
        -b data/variant_filtering/depth/${pop}.rd_filter.bed \
    > ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.rd_fil.vcf

    echo "removing non-variant SNPs from lpa-${pop} population pair"
    bcftools view -e "INFO/AF=1.00 | INFO/AF=0.00" \
        ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.rd_fil.vcf \
    > ${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf
done
```