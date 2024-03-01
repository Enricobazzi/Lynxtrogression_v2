#!/bin/bash
#SBATCH --time=1-12:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=20

# Remember to prepare the reference genome in advance! Like this:
# ${bwa} index ${reference_genome}
# ${samtools} faidx ${reference_genome}
# ${samtools} dict ${reference_genome} -o ${reference_genome/.fa/.dict}

# prepare the environment in CESGA ft3
module load cesga/2020 gcccore/system 
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.25.5
module load gatk/3.7-0-gcfedb67

# print start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Start Time: $start_time"

# the script takes as first positional argument a yaml file:
yaml=${1}

# the yaml file should have the following entries (with an example given for their values)

# sample_name: sample1
# fastq_id: 'fastq_1 fastq_2'
# fastq_r1: 'sample1_fastq_1_r1.fq sample1_fastq_2_r1.fq'
# fastq_r2: 'sample1_fastq_1_r2.fq sample1_fastq_2_r2.fq'
# in_fastq_folder: /path/to/fastqs
# reference_genome: /path/to/ref_genome
# out_bam_folder: /path/to/bams
# threads: 20
# alignment_name: align_name
# bwa: bwa
# samtools: samtools
# picard: 'java -jar ${EBROOTPICARD}/picard.jar'
# gatk: 'java -jar ${EBROOTGATK}/GenomeAnalysisTK.jar'

# function to parse the yaml file
parse_yaml() {
    local yaml_file="$1"
    local grep_word="$2"

    if [[ ! -f "$yaml_file" ]]; then
        echo "Error: YAML file '$yaml_file' not found."
        return 1
    fi

    local grep_result=$(grep "$grep_word" "$yaml_file")

    if [[ -z "$grep_result" ]]; then
        echo "Error: Word '$grep_word' not found in YAML file."
        return 1
    fi

    local trimmed_result=$(echo "$grep_result" | cut -d':' -f2 | sed -e 's/^[[:space:]]*//' | sed -e "s/^'//" -e "s/'$//")

    echo "$trimmed_result"
}

# read variables from yaml file
sample_name=$(parse_yaml $yaml "sample_name")
in_fastq_folder=$(parse_yaml $yaml "in_fastq_folder")
fastq_id=($(parse_yaml $yaml "fastq_id"))
fastq_r1=($(parse_yaml $yaml "fastq_r1"))
fastq_r2=($(parse_yaml $yaml "fastq_r2"))
reference_genome=$(parse_yaml $yaml "reference_genome")
out_bam_folder=$(parse_yaml $yaml "out_bam_folder")
threads=$(parse_yaml $yaml "threads")
alignment_name=$(parse_yaml $yaml "alignment_name")
bwa=$(parse_yaml $yaml "bwa")
samtools=$(parse_yaml $yaml "samtools")
picard=$(parse_yaml $yaml "picard")
gatk=$(parse_yaml $yaml "gatk")

for (( i = 0; i < ${#fastq_id[@]}; i++ )); do
    # Alignment
    echo " - Aligning ${fastq_id[i]} of ${sample_name} -"
    ${bwa} mem ${reference_genome} \
        ${in_fastq_folder}/${fastq_r1[i]} \
        ${in_fastq_folder}/${fastq_r2[i]} \
        -t ${threads} | 
      ${samtools} view -hbS -@ ${threads} - \
        -o ${out_bam_folder}/${sample_name}_${fastq_id[i]}_${alignment_name}.bam
    
    # Sorting
    echo " - Sorting ${sample_name}_${fastq_id[i]}_${alignment_name}.bam -"
    ${samtools} sort \
        -@ ${threads} \
        ${out_bam_folder}/${sample_name}_${fastq_id[i]}_${alignment_name}.bam \
        -o ${out_bam_folder}/${sample_name}_${fastq_id[i]}_${alignment_name}.sorted.bam
    
    rm ${out_bam_folder}/${sample_name}_${fastq_id[i]}_${alignment_name}.bam
    
    # Add Read Groups
    echo " - Adding Read Groups of ${sample_name}_${fastq_id[i]}_${alignment_name} -"
    ${picard} AddOrReplaceReadGroups \
        I=${out_bam_folder}/${sample_name}_${fastq_id[i]}_${alignment_name}.sorted.bam \
        O=${out_bam_folder}/${sample_name}_${fastq_id[i]}_${alignment_name}.sorted.rg.bam \
        RGID=${fastq_id[i]} RGLB=${sample_name}_lib \
        RGPL=Illumina RGPU=${fastq_id[i]} RGSM=${sample_name} \
        VALIDATION_STRINGENCY=SILENT
    
    rm ${out_bam_folder}/${sample_name}_${fastq_id[i]}_${alignment_name}.sorted.bam
done

if [ "${#fastq_id[@]}" -gt 1 ]; then
    
    echo "${sample_name} has multiple r1-r2 fastq pairs which need to be merged"
    
    # create list of bams to merge
    ls ${out_bam_folder}/${sample_name}_*_${alignment_name}.sorted.rg.bam \
         > ${out_bam_folder}/${sample_name}_${alignment_name}.bam.list
    
    # merge all bams of sample
    echo " - Merging ${sample_name}_${alignment_name} bams"
    ${samtools} merge -@ ${threads} \
        -r  -b ${out_bam_folder}/${sample_name}_${alignment_name}.bam.list \
        ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged.bam

    for bam in $(cat ${out_bam_folder}/${sample_name}_${alignment_name}.bam.list); do
        echo " - Removing ${bam} -"
        rm ${bam}
    done
    
    # sort the merged bam
    echo " - Sorting ${sample_name}_${alignment_name} merged bam -"
    ${samtools} sort  -@ ${threads} ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged.bam \
        -o ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.bam \
    
    rm ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged.bam

else

    echo "${sample_name} only has one r1-r2 fastq pair with id ${fastq_id[0]}"
    # Add merged_sorted to name to match alignments with multiple fastq r1-r2 pairs
    mv ${out_bam_folder}/${sample_name}_${fastq_id[0]}_${alignment_name}.sorted.rg.bam \
        ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.bam

fi


# Mark Duplicated Reads
echo " - Mark Duplicated Reads of ${sample_name}_${alignment_name} -"
${picard} MarkDuplicates \
    METRICS_FILE=${out_bam_folder}/${sample_name}_${alignment_name}.rmdup.txt \
    I=${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.bam \
    O=${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.bam \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800

rm ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.bam

# Indexing for GATK
echo " - Indexing ${sample_name}_${alignment_name} -"
${samtools} index ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.bam

# Create Target for Realignment:
echo " - RealignerTargetCreator on ${sample_name}_${fastq_id}_${alignment_name} -"
${gatk} -T RealignerTargetCreator \
    -nt ${threads} -R ${reference_genome} \
    -I ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.bam \
    -o ${out_bam_folder}/${sample_name}_realignertargetcreator.intervals

# Realign INDELs 
echo " - IndelRealigner of ${sample_name}_${alignment_name} -"
${gatk} -T IndelRealigner \
    -R ${reference_genome} \
    -targetIntervals ${out_bam_folder}/${sample_name}_realignertargetcreator.intervals \
    -I ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.bam \
    -o ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.indelrealigner.bam

rm ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.bam

# Index for downstream
echo " - Indexing ${sample_name}_${alignment_name} final BAM for downstream analyses"
${samtools} index ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.indelrealigner.bam

# Remove old indexes
rm ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.indelrealigner.bai
rm ${out_bam_folder}/${sample_name}_${alignment_name}.sorted.rg.merged_sorted.rmdup.bam.bai

# print finish time
finish_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Finish Time: $finish_time"

# write start and finish times to file
echo "Start Time: $start_time" > slurm-${SLURM_JOB_ID}.time
echo "Finish Time: $finish_time" >> slurm-${SLURM_JOB_ID}.time
