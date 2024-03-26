## Variant Calling from aligned reads

I made a list of the bams to be included in the calling:
```
# bam folder:
bam_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_bams

# create a bamlist of samples from populations we want to include in our analysis
ls ${bam_dir}/*er.bam | grep -E "ll|lp" | grep -E "ca|sm|ya|vl|ki|ur" | 
  grep -vE "ca_0249|ca_0253|sm_0138|sm_0140|sm_0185|sm_0186|sm_0221|sm_0298|sm_0359" \
  > data/bamlists/lp_ll_introgression.bamlist
```

To generate a VCF file from these BAMs, we performed variant calling using [GATK v4.2.6.1](https://gatk.broadinstitute.org/hc/en-us), using the -ERC GVCF option, to generate genome VCFs of each sample. The [call_gvcf_ref_bam_outgvcf_bed](src/calling/call_gvcf_ref_bam_outgvcf_bed.sh) script performs the calling on a single chromosome (or any region provided as a BED file), so that jobs can be faster and parallelized. Each chromosome selected is further divided in 8 parts, so I can run simultaneously the 8 commands in one node.

This is how I selected the chromosomes to be used and divided them in 8 using [bedtools makewindows](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/bedtools-makewindows):
```
# define reference genome folder
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
# create a folder to store each chromosome's bed file
mkdir ${ref_dir}/chromosome_beds
# create a list of autosomic chromosomes
grep -i "chr" ${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa.fai | grep -viE "chry|chrx" | cut -f1 > ${ref_dir}/autosomic_scaffolds_list.txt
# create each chromosome's bed in a loop
for chr in $(cat ${ref_dir}/autosomic_scaffolds_list.txt); do
    grep ${chr} ${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa.fai |
        awk '{print $1, 0, $2}' | tr ' ' '\t' |
        bedtools makewindows -b - -n 8 \
        > ${ref_dir}/chromosome_beds/${chr}.mk8win.bed
    for n in {1..8}; do
        sed "${n}q;d" ${ref_dir}/chromosome_beds/${chr}.mk8win.bed \
            > ${ref_dir}/chromosome_beds/${chr}.${n}.bed
    done
done
```

Calling of each chromosome of each sample was performed using the [sbatch_haplotypecaller_ref_inbam_gvcfdir_chr](src/calling/sbatch_haplotypecaller_ref_inbam_gvcfdir_chr.sh) script, where variants found in the 8 sub-chromosome windows are called in parallel.

This was sbatched to the ft3 queue as follows:
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
ref=${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa
inbams=($(cat data/bamlists/lp_ll_introgression.bamlist))
chromosomes=($(cat ${ref_dir}/autosomic_scaffolds_list.txt))
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_gvcfs

# from 0 to 72!
for i in {70..72}; do
    inbam=${inbams[${i}]}
    sample=$(basename -a ${inbam} | cut -d'_' -f1,2,3,4)

    for chr in ${chromosomes[*]}; do
        echo "call $chr of $sample"
        sbatch \
            --job-name=${sample}.${chr} \
            --output=logs/calling/gvcf.${sample}.${chr}.out \
            --error=logs/calling/gvcf.${sample}.${chr}.err \
            src/calling/sbatch_haplotypecaller_ref_inbam_gvcfdir_chr.sh \
            ${ref} \
            ${inbam} \
            ${gvcf_dir} \
            ${chr}
    done
done

# find any failed runs:
for err in $(ls logs/calling/gvcf.*.err); do
    ndone=$(grep "HaplotypeCaller done" ${err} | wc -l)
    if [ "${ndone}" -lt 8 ]; then
        echo "${err} has ${ndone}"
    fi
done
```

A gvcf of each chromosome of each sample is generated using the [sbatch_mergevcfs_sample_gvcfdir_chrlist](src/calling/sbatch_mergevcfs_sample_gvcfdir_chrlist.sh) script by merging each sub-chromosome window gvcf.

This was sbatched to the ft3 queue as follows:
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
ref=${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa
samples=($(cat data/sample.list))
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_gvcfs
chr_list=${ref_dir}/autosomic_scaffolds_list.txt

for sample in ${samples[*]}; do
    echo "sbatch sbatch_mergevcfs_sample_gvcfdir_chrlist of $sample"
    sbatch \
        --job-name=${sample}.mergegvcfs \
        --output=logs/calling/mergegvcfs.${sample}.out \
        --error=logs/calling/mergegvcfs.${sample}.err \
        src/calling/sbatch_mergevcfs_sample_gvcfdir_chrlist.sh \
        ${sample} \
        ${gvcf_dir} \
        ${chr_list}
done

# to eliminate the sub-chromosome windows gvcfs when completed:
for sample in ${samples[*]}; do for chr in $(cat $chr_list); do rm ${gvcf_dir}/${sample}.${chr}.*.g.vcf; done; done
```

Then a gvcf for each chromosome containing all samples of the project is created for each chromosome using the [sbatch_combinegvcfs_gvcflist_ref_outgvcf](src/calling/sbatch_combinegvcfs_gvcflist_ref_outgvcf.sh) script.

This was sbatched to the ft3 queue as follows:
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
ref=${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa
samples=($(cat data/sample.list))
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_gvcfs
chr_list=${ref_dir}/autosomic_scaffolds_list.txt

for chr in $(cat ${chr_list}); do
    ls ${gvcf_dir}/*.${chr}.g.vcf > tmp.${chr}.gvcf.list
    outgvcf=${gvcf_dir}/lynxtrogression_v2.${chr}.g.vcf
    echo "sbatch sbatch_combinegvcfs_gvcflist_ref_outgvcf ${outgvcf}"
    sbatch \
        --job-name=${chr}.combinegvcfs \
        --output=logs/calling/combinegvcfs.${chr}.out \
        --error=logs/calling/combinegvcfs.${chr}.err \
        src/calling/sbatch_combinegvcfs_gvcflist_ref_outgvcf.sh \
        tmp.${chr}.gvcf.list \
        ${ref} \
        ${outgvcf}
done

# when finished clean tmp files
rm tmp.mLynRuf2.2_Chr*.list
```

To generate vcf files from the chromosome gvcfs I run the [sbatch_genotypegvcfs_ref_ingvcf_outvcf](src/calling/sbatch_genotypegvcfs_ref_ingvcf_outvcf.sh) script.

This was sbatched to the ft3 queue as follows:
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
ref=${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa
chr_list=${ref_dir}/autosomic_scaffolds_list.txt
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_gvcfs
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs

for chr in $(cat ${chr_list}); do
    ingvcf=${gvcf_dir}/lynxtrogression_v2.${chr}.g.vcf
    outvcf=${vcf_dir}/lynxtrogression_v2.${chr}.vcf
    echo "sbatch genotypegvcf ${ingvcf}"
    sbatch \
        --job-name=${chr}.genotypegvcfs \
        --output=logs/calling/genotypegvcfs.${chr}.out \
        --error=logs/calling/genotypegvcfs.${chr}.err \
        src/calling/sbatch_genotypegvcfs_ref_ingvcf_outvcf.sh \
        ${ref} \
        ${ingvcf} \
        ${outvcf}
done
```

Finally to create a vcf with all the chromosomes of all samples I run the [sbatch_bcftoolsconcat_vcflist_outvcf](src/calling/sbatch_bcftoolsconcat_vcflist_outvcf.sh) script.

This was sbatched to the ft3 queue as follows:
```
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_vcfs
ls ${vcf_dir}/lynxtrogression_v2.*.vcf > tmp.vcf.list
vcflist=tmp.vcf.list
outvcf=${vcf_dir}/lynxtrogression_v2.autosomic_scaffolds.vcf.gz

sbatch src/calling/sbatch_bcftoolsconcat_vcflist_outvcf.sh \
    ${vcflist} ${outvcf}

rm tmp.vcf.list
```
