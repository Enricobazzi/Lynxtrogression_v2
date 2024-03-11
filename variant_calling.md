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

Calling of each chromosome of each sample was sbatched to the ft3 queue as follows:
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_rufus_mLynRuf2.2
ref=${ref_dir}/mLynRuf2.2.revcomp.scaffolds.fa
inbams=($(cat data/bamlists/lp_ll_introgression.bamlist))
chromosomes=($(cat ${ref_dir}/autosomic_scaffolds_list.txt))
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynRuf2.2_ref_gvcfs

# to 73!
for i in {10..14}; do
    inbam=${inbams[${i}]}
    sample=$(basename -a ${inbam} | cut -d'_' -f1,2,3,4)

    for chr in ${chromosomes[*]}; do
        echo "call $chr of $sample"
        sbatch src/calling/sbatch_haplotypecaller_ref_inbam_gvcfdir_chr.sh \
            ${ref} \
            ${inbam} \
            ${gvcf_dir} \
            ${chr}
    done
done
```
