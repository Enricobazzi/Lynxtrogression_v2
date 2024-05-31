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

### Preparations

To run GADMA2 I need 3 things:
- vcf file: I will use the pruned vcf of each population pair
- popmap file: which assigns a population to each sample in the vcf
- parameter file: where the gadma run specifications are stored, I have one for each population pair

In the parameter file I selected the following based on the population pair:
```
Input data: <pop-pair.vcf>,<pop-pair.popmap>
Projections: [<nhaplo_pop1>, <nhaplo_pop2>]
Sequence length: <pop-pair_genome_length_after_filters>
```
Here pop1 and pop2 are in the order they appear in the popmap file (in my case always pop1 = Eurasian lynx and pop2 = Iberian lynx) and their projections is the number of haplotypes they have (number of individuals * 2) because while pruning I eliminated all missing data.

Other parameters are the same regardless of population pair
```
Outgroup: False
Mutation rate: 6e-9
Initial structure: [2, 2]
Final structure: [2, 3]
Number of repeats: 10
Number of processes: 10
Split fractions: False
Lower bound of first split: 20000
Dynamics: Sud, Exp
Linked SNP's: True
Min_N: 0.001
Max_N: 50
Max_T: 10
Pts: [50, 60, 70]
Size_of_generation: 10
# Fractions
N_elitism: 3
P_mutation: 0.2
P_crossover: 0.3
P_random: 0.2
Mean_mutation_strength: 0.776
Const_for_mutation_strength: 1.302
Mean_mutation_rate: 0.273
Const_for_mutation_rate: 1.475
```

### Running GADMA2 *to be completed*

Run gadma:
```
pair=lpa-wel
for n in {1..50}; do
    sbatch --job-name=${n}_gadma_${pair} \
        --output=logs/demographic_inference/${pair}.gadma.${n}.out \
        --error=logs/demographic_inference/${pair}.gadma.${n}.err \
        src/demographic_inference/run_gadma.sh ${pair} ${n} ${pair}_${n}
done
```

Resume the unfinished runs:
```
# first resumed
pair=lpa-wel
for n in {1..50}; do
    ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}/GADMA.log | wc -l)
    if [ $ncomplete -ne 10 ]; then
        echo "RESUMING RUN ${pair}_${n}:"
        sbatch --job-name=${n}_gadma_${pair} \
            --output=logs/demographic_inference/${pair}.gadma.${n}.out \
            --error=logs/demographic_inference/${pair}.gadma.${n}.err \
            src/demographic_inference/run_gadma.sh ${pair} ${n} ${pair}_${n}
    fi
done

pair=lpa-wel
for n in {1..50}; do
    if [ -d data/demographic_inference/${pair}_${n}_resumed ]; then
    ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}_resumed/GADMA.log | wc -l)
        if [ $ncomplete -eq 10 ]; then
            echo "RUN ${pair}_${n}_resumed COMPLETE!"
        fi
    fi
done

# second resumed
pair=lpa-wel
for n in {1..50}; do
    ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}_resumed/GADMA.log | wc -l)
    if [ $ncomplete -ne 10 ]; then
        echo "RESUMING RUN ${pair}_${n}_resumed :"
        sbatch --job-name=${n}_gadma_${pair} \
            --output=logs/demographic_inference/${pair}.gadma.${n}.out \
            --error=logs/demographic_inference/${pair}.gadma.${n}.err \
            src/demographic_inference/run_gadma.sh ${pair} ${n} ${pair}_${n}_resumed
    fi
done

pair=lpa-wel
for n in {1..50}; do
    if [ -d data/demographic_inference/${pair}_${n}_resumed_resumed ]; then
    ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}_resumed_resumed/GADMA.log | wc -l)
        if [ $ncomplete -eq 10 ]; then
            echo "RUN ${pair}_${n}_resumed_resumed COMPLETE!"
        else
            echo "RUN ${pair}_${n}_resumed_resumed: ${ncomplete} completed"
        fi
    fi
done

# third resumed
pair=lpa-wel
for n in {1..50}; do
    ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}_resumed_resumed/GADMA.log | wc -l)
    if [ $ncomplete -ne 10 ]; then
        echo "RESUMING RUN ${pair}_${n}_resumed_resumed :"
        sbatch --job-name=${n}_gadma_${pair} \
            --output=logs/demographic_inference/${pair}.gadma.${n}.out \
            --error=logs/demographic_inference/${pair}.gadma.${n}.err \
            src/demographic_inference/run_gadma.sh ${pair} ${n} ${pair}_${n}_resumed_resumed
    fi
done

pair=lpa-wel
for n in {1..50}; do
    if [ -d data/demographic_inference/${pair}_${n}_resumed_resumed_resumed ]; then
    ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}_resumed_resumed_resumed/GADMA.log | wc -l)
        if [ $ncomplete -eq 10 ]; then
            echo "RUN ${pair}_${n}_resumed_resumed_resumed COMPLETE!"
        else
            echo "RUN ${pair}_${n}_resumed_resumed_resumed: ${ncomplete} completed"
        fi
    fi
done

# fourth resumed (gave it 3 days)
pair=lpa-wel
for n in {1..50}; do
    if [ -d data/demographic_inference/${pair}_${n}_resumed_resumed_resumed ]; then
        ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}_resumed_resumed_resumed/GADMA.log | wc -l)
        if [ $ncomplete -ne 10 ]; then
            echo "RESUMING RUN ${pair}_${n}_resumed_resumed_resumed :"
            sbatch --job-name=${n}_gadma_${pair} \
                --output=logs/demographic_inference/${pair}.gadma.${n}.out \
                --error=logs/demographic_inference/${pair}.gadma.${n}.err \
                src/demographic_inference/run_gadma.sh ${pair} ${n} ${pair}_${n}_resumed_resumed_resumed
        fi
    fi
done

pair=lpa-wel
for n in {1..50}; do
    if [ -d data/demographic_inference/${pair}_${n}_resumed_resumed_resumed_resumed ]; then
    ncomplete=$(grep "algor" data/demographic_inference/${pair}_${n}_resumed_resumed_resumed_resumed/GADMA.log | wc -l)
        if [ $ncomplete -eq 10 ]; then
            echo "RUN ${pair}_${n}_resumed_resumed_resumed_resumed COMPLETE!"
        else
            echo "RUN ${pair}_${n}_resumed_resumed_resumed_resumed: ${ncomplete} completed"
        fi
    fi
done

```

Get results (demes yamls and moments scripts) *explore all the resumed folders*
```
pair=lpa-wel

# DEMES YAMLS = final_best_logLL_model_demes_code.py.yml
for n in {1..50}; do
    for k in {1..10}; do
        yaml=data/demographic_inference/${pair}_${n}_resumed_resumed_resumed/${k}/final_best_logLL_model_demes_code.py.yml
        if [ -f $yaml ]; then
            cp $yaml data/demographic_inference/${pair}_best_yamls/${pair}_${n}_${k}_final_best_model.yaml
        fi
    done
done
# MOMENTS SCRIPTS = final_best_logLL_model_moments_code.py
for n in {1..50}; do
    for k in {1..10}; do
        moms=data/demographic_inference/${pair}_${n}_resumed_resumed_resumed/${k}/final_best_logLL_model_moments_code.py
        if [ -f $moms ]; then
            cp $moms data/demographic_inference/${pair}_best_moments/${pair}_${n}_${k}_final_best_moments.py
        fi
    done
done
```

```
pair=lpa-wel
# LIKELIHOOD TABLE
echo "run log-like" > ll_table.txt

for n in {1..50}; do
    for k in {1..10}; do
        moms=data/demographic_inference/${pair}_${n}/${k}/final_best_logLL_model_moments_code.py
        if [ -f $moms ]; then
            ll=$(tail data/demographic_inference/${pair}_${n}/${k}/GADMA_GA.log | grep "y:" | cut -d':' -f2 | cut -d' ' -f2 | tail -1)
            echo "${pair}_${n}_${k}: ${ll}"
            echo "${pair}_${n}_${k} ${ll}" >> ll_table.txt
        fi
    done
done

for n in {1..50}; do
    for k in {1..10}; do
        moms=data/demographic_inference/${pair}_${n}_resumed/${k}/final_best_logLL_model_moments_code.py
        if [ -f $moms ]; then
            ll=$(tail data/demographic_inference/${pair}_${n}_resumed/${k}/GADMA_GA.log | grep "y:" | cut -d':' -f2 | cut -d' ' -f2 | tail -1)
            echo "${pair}_${n}_${k}: ${ll}"
            echo "${pair}_${n}_${k} ${ll}" >> ll_table.txt
        fi
    done
done

for n in {1..50}; do
    for k in {1..10}; do
        moms=data/demographic_inference/${pair}_${n}_resumed_resumed/${k}/final_best_logLL_model_moments_code.py
        if [ -f $moms ]; then
            ll=$(tail data/demographic_inference/${pair}_${n}_resumed_resumed/${k}/GADMA_GA.log | grep "y:" | cut -d':' -f2 | cut -d' ' -f2 | tail -1)
            echo "${pair}_${n}_${k}: ${ll}"
            echo "${pair}_${n}_${k} ${ll}" >> ll_table.txt
        fi
    done
done


for n in {1..50}; do
    for k in {1..10}; do
        moms=data/demographic_inference/${pair}_${n}_resumed_resumed_resumed/${k}/final_best_logLL_model_moments_code.py
        if [ -f $moms ]; then
            ll=$(tail data/demographic_inference/${pair}_${n}_resumed_resumed_resumed/${k}/GADMA_GA.log | grep "y:" | cut -d':' -f2 | cut -d' ' -f2 | tail -1)
            echo "${pair}_${n}_${k}: ${ll}"
            echo "${pair}_${n}_${k} ${ll}" >> ll_table.txt
        fi
    done
done

sort -rnk 2 ll_table.txt | uniq > tmp && mv tmp ll_table.txt
```

### Obtaining confidence intervals from best runs *to be completed*

Instructions in the [gadma manual page](https://gadma.readthedocs.io/en/latest/user_manual/confidence_intervals.html)

In order to obtain confidence intervals for our results we need to bootstrap our data. The [make_bootstrap_from_gadma](src/demographic_inference/make_bootstrap_from_gadma.py) python script will do it using the information from the gadma parameter file. For each population pair I create 100 bootstrap replicates and save them in a folder called `<pair>_bootstrap`
```
module load cesga/system miniconda3/22.11.1-1
source activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/gadma2

for pair in lpa-wel lpa-eel lpa-sel; do
    echo "bootstrapping ${pair}"
    python src/demographic_inference/make_bootstrap_from_gadma.py \
        data/demographic_inference/${pair}.gadma.parameters \
        data/demographic_inference/${pair}_bootstrap \
        100
done
```

To run local search from initial and bootstrapped Frequency Spectrums I use the `gadma-run_ls_on_boot_data` gadma script in the [sbatch_ls_on_boot.sh](src/demographic_inference/sbatch_ls_on_boot.sh) custom script
```
pair=lpa-wel
for run in 12_9 6_2 20_7; do
    sbatch --job-name=${run}_CI \
        --output=logs/demographic_inference/${pair}_${run}.CI.out \
        --error=logs/demographic_inference/${pair}_${run}.CI.err \
        src/demographic_inference/sbatch_ls_on_boot.sh \
            data/demographic_inference/${pair}_bootstrap \
            data/demographic_inference/${pair}_best_moments/${pair}_${run}_final_best_moments.py \
            data/demographic_inference/${pair}_CI/${pair}_${run}
done
```
