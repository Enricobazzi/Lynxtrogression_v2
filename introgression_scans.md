# Introgression Scans

## Download Intronets

I need [introNets](https://github.com/SchriderLab/introNets) version of ms (msmodified) to run the simulations. I download the folder from github and compile msmodified as suggested by Dylan.
```
cd src
git clone https://github.com/SchriderLab/introNets.git
cd introNets/msmodified
gcc -o ms ms.c streec.c rand1.c -lm
```
To setup a conda environment to run its scripts:
```
## create environment
# on genomics cluster
conda create --prefix=/home/ebazzicalupo/introNets/intronets python=3.9
conda activate /home/ebazzicalupo/introNets/intronets
# on cesga cluster
conda create --prefix=/mnt/netapp1/Store_CSIC/home/csic/eye/eba/intronets python=3.9
conda activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/intronets

# install stuff:
conda install -c conda-forge mpi4py openmpi
pip install "numpy<1.25"
pip install seriate
pip install scipy
pip install scikit-learn
pip install h5py
python -m pip install -U matplotlib
pip uninstall protobuf
pip install "protobuf<3.20"
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install pandas
pip install seaborn
pip install prettytable
```

## Simulate training data

To simulate the data for training the [simulate_data.py](src/introgression_scans/simulate_data.py) takes a YAML file (`--demes_yaml`) with the demographic model in demes format (generated by GADMA2), a CSV (`--confint`) with the confidence intervals around the parameters (generated by [plot_boot_params.py](src/demographic_inference/plot_boot_params.py)), and a migration scenario (`--migration`, see below), and outputs a specified number of simulations (`--nreps`) to an output directory (`--odir`), using a msmodified simulator (`--path_to_msmodified`).

The migration scenarios to be simulated are:
- ab: migration from population 1 (eurasian lynx) to population 2 (iberian lynx) in forward time (`-es Tmig 2 Pmig -ej Tmig 3 1`)
- ba: migration from population 2 (iberian lynx) to population 1 (eurasian lynx) in forward time (`-es Tmig 1 Pmig -ej Tmig 3 2`)
- abba: bidirectional migration, with first ab and then ba in forward time
- baab: bidiractional migration, with first ba and then ab in forward time
- none: no migration

The time of migration is set to be a random time between the present and 10 thousand generations ago and the amount of migration is a random value between 0.01 and 0.3.

Both abba and baab will go into the 'bi' category for the discriminator, but are divided here because of how you need to specify the order of migrations in ms. This means we will simulate X ab, ba and none, and X/2 abba and baab for each demographic model.

```
pop_pair=lpa-wel
pop_pair=lpa-eel
pop_pair=lpa-sel

for migration in ab ba bi none; do
    mkdir data/introgression_scans/simulations/${pop_pair}_${migration}_sims/
done

for model in 12_9 20_7 6_2; do
for model in 34_7 38_4 30_1; do
for model in 12_6 18_7 18_10; do
    mkdir data/introgression_scans/simulations/${pop_pair}_ab_sims/${model}/
    mkdir data/introgression_scans/simulations/${pop_pair}_ba_sims/${model}/
    mkdir data/introgression_scans/simulations/${pop_pair}_bi_sims/${model}_abba/
    mkdir data/introgression_scans/simulations/${pop_pair}_bi_sims/${model}_baab/
    mkdir data/introgression_scans/simulations/${pop_pair}_none_sims/${model}/
done

for model in 12_9 20_7 6_2; do
for model in 34_7 38_4 30_1; do
for model in 12_6 18_7 18_10; do
    for migration in ab ba none; do
        python src/introgression_scans/simulate_data.py \
            --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
            --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
            --path_to_msmodified src/introNets/msmodified/ms \
            --migration ${migration} \
            --nreps 20000 \
            --odir data/introgression_scans/simulations/${pop_pair}_${migration}_sims/${model}/
    done
done

pop_pair=lpa-sel
migration=none
conda activate lynxtrogression_v2

for model in 12_6 18_7 18_10; do
python src/introgression_scans/simulate_data.py \
    --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
    --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
    --path_to_msmodified src/introNets/msmodified/ms \
    --migration ${migration} \
    --nreps 20000 \
    --odir data/introgression_scans/simulations/${pop_pair}_${migration}_sims/${model}/
done

for model in 12_9 20_7 6_2; do
for model in 34_7 38_4 30_1; do
for model in 12_6 18_7 18_10; do
    for migration in abba baab; do
        python src/introgression_scans/simulate_data.py \
            --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
            --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
            --path_to_msmodified src/introNets/msmodified/ms \
            --migration ${migration} \
            --nreps 10000 \
            --odir data/introgression_scans/simulations/${pop_pair}_bi_sims/${model}_${migration}/
    done
done
```

## Filter simulated data for training

I filter out from simulations the ones with < 128 segsites using [filter_sims.py](src/introgression_scans/filter_sims.py), or [introNets format.py](https://github.com/SchriderLab/introNets/blob/main/src/data/format.py) will throw errors and corrupt the hdf5 file:
```
pop_pair=lpa-wel
pop_pair=lpa-eel
pop_pair=lpa-sel

for migration in ab ba bi none; do
    mkdir data/introgression_scans/simulations/${pop_pair}_${migration}_filtered_sims/
done

for model in 12_9 20_7 6_2; do
for model in 34_7 38_4 30_1; do
for model in 12_6 18_7 18_10; do
    mkdir data/introgression_scans/simulations/${pop_pair}_ab_filtered_sims/${model}/
    mkdir data/introgression_scans/simulations/${pop_pair}_ba_filtered_sims/${model}/
    mkdir data/introgression_scans/simulations/${pop_pair}_bi_filtered_sims/${model}_abba/
    mkdir data/introgression_scans/simulations/${pop_pair}_bi_filtered_sims/${model}_baab/
    mkdir data/introgression_scans/simulations/${pop_pair}_none_filtered_sims/${model}/
done

for model in 12_9 20_7 6_2; do
for model in 34_7 38_4 30_1; do
for model in 12_6 18_7 18_10; do
    for migration in ab ba none; do
        echo "filtering ${pop_pair}_${migration} of ${model}"
        python src/introgression_scans/filter_sims.py \
            --idir data/introgression_scans/simulations/${pop_pair}_${migration}_sims/${model}/ \
            --odir data/introgression_scans/simulations/${pop_pair}_${migration}_filtered_sims/${model}/ \
            --n_sites 128 --n_sims 12000
    done
done

for model in 12_9 20_7 6_2; do
for model in 34_7 38_4 30_1; do
for model in 12_6 18_7 18_10; do
    for migration in abba baab; do
        echo "filtering ${pop_pair}_${migration} of ${model}"
        python src/introgression_scans/filter_sims.py \
            --idir data/introgression_scans/simulations/${pop_pair}_bi_sims/${model}_${migration}/ \
            --odir data/introgression_scans/simulations/${pop_pair}_bi_filtered_sims/${model}_${migration}/ \
            --n_sites 128 --n_sims 6000
    done
done
```

## Format simulations

To format my simulations for training I run [introNets format.py](https://github.com/SchriderLab/introNets/blob/main/src/data/format.py):
```
## ADD HDF5 FOLDER TO OFOLDER

conda activate ~/introNets/intronets
pop_pair=lpa-wel
pop_sizes="40,44"
conda activate ~/introNets/intronets
pop_pair=lpa-eel
pop_sizes="38,44"
conda activate ~/introNets/intronets
pop_pair=lpa-sel
pop_sizes="24,44"

mkdir data/introgression_scans/simulations/${pop_pair}_hdf5s

conda activate ~/introNets/intronets
pop_pair=lpa-sel
pop_sizes="24,44"
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/introgression_scans/simulations/${pop_pair}_ab_filtered_sims/ \
    --ofile data/introgression_scans/simulations/${pop_pair}_hdf5s/${pop_pair}_ab.hdf5 \
    --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 1 |& tee logs/introgression_scans/format_${pop_pair}_ab.log

conda activate ~/introNets/intronets
pop_pair=lpa-sel
pop_sizes="24,44"
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/introgression_scans/simulations/${pop_pair}_ba_filtered_sims/ \
    --ofile data/introgression_scans/simulations/${pop_pair}_hdf5s/${pop_pair}_ba.hdf5 \
    --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 0 |& tee logs/introgression_scans/format_${pop_pair}_ba.log

conda activate ~/introNets/intronets
pop_pair=lpa-sel
pop_sizes="24,44"
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/introgression_scans/simulations/${pop_pair}_bi_filtered_sims/ \
    --ofile data/introgression_scans/simulations/${pop_pair}_hdf5s/${pop_pair}_bi.hdf5 \
    --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop -1 |& tee logs/introgression_scans/format_${pop_pair}_bi.log

conda activate ~/introNets/intronets
pop_pair=lpa-sel
pop_sizes="24,44"
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/introgression_scans/simulations/${pop_pair}_none_filtered_sims/ \
    --ofile data/introgression_scans/simulations/${pop_pair}_hdf5s/${pop_pair}_none.hdf5 \
    --pop_sizes ${pop_sizes} --out_shape 2,44,128 --include_zeros |& tee logs/introgression_scans/format_${pop_pair}_none.log
```

## Check formatted simulations

I can create a fasta file to visually check the formatted simulations to see if they make sense using [write_fastas_from_hdf5.py](src/introgression_scans/write_fastas_from_hdf5.py):
```
conda activate ~/introNets/intronets
pop_pair=lpa-wel
pop_pair=lpa-eel

# mkdir data/introgression_scans/simulations/fastas/
for mig in ab ba bi none; do
    mkdir data/introgression_scans/simulations/fastas/${pop_pair}_${mig}
done

for mig in ab ba bi none; do
    python src/introgression_scans/write_fastas_from_hdf5.py \
        --ifile data/introgression_scans/simulations/${pop_pair}_hdf5s/${pop_pair}_${mig}.hdf5 \
        --odir data/introgression_scans/simulations/fastas/${pop_pair}_${mig} \
        --nseqs 25
done
```

## Train a discriminator model

To train a discriminator of my 4 classes of introgression I run [introNets train_discriminator.py](https://github.com/SchriderLab/introNets/blob/main/src/models/train_discriminator.py):

```
pop_pair=lpa-wel
script logs/introgression_scans/${pop_pair}_train_disc.log
conda activate ~/introNets/intronets
pop_pair=lpa-wel
taskset -c 1,2,3,4,5,6,7,8,9,10 \
    python src/introNets/src/models/train_discriminator.py \
        --idir data/introgression_scans/simulations/${pop_pair}_hdf5s/ --odir data/introgression_scans/${pop_pair}_discriminator/ --n_classes 4

pop_pair=lpa-eel
script logs/introgression_scans/${pop_pair}_train_disc.log
conda activate ~/introNets/intronets
pop_pair=lpa-eel
taskset -c 11,12,13,14,15,16,17,18,19,20 \
    python src/introNets/src/models/train_discriminator.py \
        --idir data/introgression_scans/simulations/${pop_pair}_hdf5s/ --odir data/introgression_scans/${pop_pair}_discriminator/ --n_classes 4

pop_pair=lpa-sel
script logs/introgression_scans/${pop_pair}_train_disc.log
conda activate ~/introNets/intronets
pop_pair=lpa-sel
mkdir data/introgression_scans/${pop_pair}_discriminator/
taskset -c 1,2,3,4,5,6,7,8,9,10 \
    python src/introNets/src/models/train_discriminator.py \
        --idir data/introgression_scans/simulations/${pop_pair}_hdf5s/ --odir data/introgression_scans/${pop_pair}_discriminator/ --n_classes 4

```

## Evaluate the discriminator

TBD

## Real data preparation

To transform the real data (phased vcf) to a hdf5 file that can be analyzed by the discriminator I first convert it to numpy's [NPZ](https://imageio.readthedocs.io/en/v2.5.0/format_npz.html) format. For this I use the [phasedVcfToNpz.py](src/introgression_scans/phasedVcfToNpz.py) script wrote by Dan. When using this script the `species1` will be saved as `sechHeader` + `sechMatrix` and `species2` as `simHeader` + `simMatrix` in the NPZ (important for hdf5 formatting below). One separate NPZ file is generated for each chromosome:
```
# vcfdir
vcfdir=/GRUPOS/grupolince/mLynRuf2.2_ref_vcfs
# refdir
refdir=/GRUPOS/grupolince/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))
# eurasian pop
pop=wel
pop=eel
pop=sel
# vcf
vcfFileName=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf
# species1 - eurasian lynx
species1ListFileName=data/${pop}.list
# species2 - iberian lynx
species2ListFileName=data/lpa.list
# reference genome (masked if needed - not our case)
maskedRefFileName=${refdir}/mLynRuf2.2.revcomp.scaffolds.fa

for chr in ${chrs[@]}; do
    # output file
    npzFileName=data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz
    echo "generating ${npzFileName}"
    # run script
    python src/introgression_scans/phasedVcfToNpz.py \
        $vcfFileName $species1ListFileName $species2ListFileName $maskedRefFileName $chr $npzFileName
done
```

I use [introNets format_npz.py](https://github.com/SchriderLab/introNets/blob/main/src/data/format_npz.py) to convert these files to hdf5. I set `--keys` to `sechMatrix,simMatrix` to have the eurasian lynx population as population1 and the iberian lynx population as population2.
```
conda activate ~/introNets/intronets

# refdir
refdir=/GRUPOS/grupolince/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))
# eurasian pop
pop=wel

for chr in ${chrs[@]}; do
    in_npz=data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz
    out_hdf5=data/introgression_scans/hdf5_files/lpa-${pop}.${chr}.hdf5
    echo "generating ${out_hdf5}"
    mpirun -n 10 python src/introNets/src/data/format_npz.py \
        --verbose \
        --ifile ${in_npz} \
        --ofile ${out_hdf5} \
        --pop_sizes 40,44 --out_shape 2,44,128 \
        --keys sechMatrix,simMatrix
done
```

## NOT CURRENTLY WORKING - Apply discriminator to real data HDF5 - 

*THIS IS NOT WORKING AS INTENDED! THERE IS A BUG IN FORMAT_NPZ.PY*

To apply the discriminator model to the real data I wrote a modified version of [introNets apply_disc.py](https://github.com/SchriderLab/introNets/blob/main/src/models/apply_disc.py) which I located in the introNets folder for easier access to required imports. I put a copy of my version ([apply_disc_eb.py](src/introgression_scans/apply_disc_eb.py)) in this repo, but in order to work properly it has to be located in the `src/models/` folder of introNets.

```
conda activate ~/introNets/intronets

ifile=data/introgression_scans/hdf5_files/lpa-wel.mLynRuf2.2_ChrA1.hdf5
ofile=data/introgression_scans/lpa-wel_predictions/mLynRuf2.2_ChrA1.predictions.new.csv

python src/introNets/src/models/apply_disc_eb.py \
    --weights data/introgression_scans/discriminator/test.weights \
    --ifile ${ifile} \
    --ofile ${ofile}
```

## Apply discriminator to real data NPZ

Since there is a problem when generating the hdf5 file from npz data, until it's fixed I will use different script, [apply_disc_to_npz.py](src/introgression_scans/apply_disc_to_npz.py), that applies the model directly to the NPZ. This is not ideal since it's much slower in sorting the haplotypes in the windows because it's not being parallelized as in `format_npz.py`. I use it until the other is fixed. I placed a copy of it in the `src/models/` folder of introNets so it can load the model properly.

```
conda activate ~/introNets/intronets
pop_pair=lpa-wel
pop_sizes="40,44"
conda activate ~/introNets/intronets
pop_pair=lpa-eel
pop_sizes="38,44"
conda activate ~/introNets/intronets
pop_pair=lpa-sel
pop_sizes="24,44"

mkdir data/introgression_scans/lpa-${pop}_predictions/

chr=mLynRuf2.2_ChrA1
chr=mLynRuf2.2_ChrC1
chr=mLynRuf2.2_ChrB1
chr=mLynRuf2.2_ChrA2_rc
chr=mLynRuf2.2_ChrC2
chr=mLynRuf2.2_ChrB2_rc
chr=mLynRuf2.2_ChrB3
chr=mLynRuf2.2_ChrB4_rc
chr=mLynRuf2.2_ChrA3_rc
chr=mLynRuf2.2_ChrD1
chr=mLynRuf2.2_ChrD4
chr=mLynRuf2.2_ChrD3
chr=mLynRuf2.2_ChrD2
chr=mLynRuf2.2_ChrF2
chr=mLynRuf2.2_ChrF1_rc
chr=mLynRuf2.2_ChrE2_rc
chr=mLynRuf2.2_ChrE1
chr=mLynRuf2.2_ChrE3_rc

python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4
```

## Create fastas from VCF to check predictions

To manually check predicted regions with introgression I run [src/introgression_scans/write_fastas_from_vcf.py](src/introgression_scans/write_fastas_from_vcf.py), then I can check them with an alignment viewer (I use [Jalview](https://www.jalview.org/)):
```
# refdir
refdir=/GRUPOS/grupolince/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))
# vcfdir
vcfdir=/GRUPOS/grupolince/mLynRuf2.2_ref_vcfs
# eurasian pop
pop=wel
pop=eel
pop=sel
# ivcf
ivcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf

mkdir data/introgression_scans/lpa-${pop}_vcf_fastas/

for chr in ${chrs[@]}; do
    echo "fastas of ${chr}"
    mkdir data/introgression_scans/lpa-${pop}_vcf_fastas/${chr}
    python src/introgression_scans/write_fastas_from_vcf.py \
        --ivcf ${ivcf} \
        --chr ${chr} \
        --odir data/introgression_scans/lpa-${pop}_vcf_fastas/${chr} \
        --win_size 128 \
        --step_size 64
done
```

From NPZ

```
conda activate ~/introNets/intronets
pop=wel
pop_sizes="40,44"
conda activate ~/introNets/intronets
pop=eel
pop_sizes="38,44"
conda activate ~/introNets/intronets
pop=sel
pop_sizes="24,44"

# refdir
refdir=/GRUPOS/grupolince/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))

mkdir data/introgression_scans/lpa-${pop}_binary_fastas/

for chr in ${chrs[@]}; do
    echo "writing binary fastas of ${chr}"
    # mkdir data/introgression_scans/lpa-${pop}_binary_fastas/${chr}
    python src/introNets/src/models/write_binary_fastas_from_npz.py \
        --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
        --odir data/introgression_scans/lpa-${pop}_binary_fastas/${chr} \
        --pop_sizes ${pop_sizes} \
        --shape 2,44,128 \
        --step_size 64
done
```

## Transform the probabilities to predictions and save them to AB, BA, and BI introgressed regions bed files *to be finished*

Create beds
```
# vcfdir
vcfdir=/GRUPOS/grupolince/mLynRuf2.2_ref_vcfs
# eurasian pop
pop=wel
pop=eel
pop=sel
# ivcf
ivcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf

python src/introgression_scans/get_beds_from_preds.py \
    --ivcf ${ivcf} \
    --folder data/introgression_scans/lpa-${pop}_predictions/ \
    --pthresh 0.95 \
    --wsize 128 \
    --step 64
```
Merge consecutive windows
```
pop=wel
pop=eel
pop=sel

ls data/introgression_scans/lpa-${pop}_predictions/*095.bed

bedtools merge \
    -i data/introgression_scans/lpa-${pop}_predictions/ab_introgressed_095.bed \
    > data/introgression_scans/lpa-${pop}_predictions/ab_introgressed_095.merged.bed

bedtools merge \
    -i data/introgression_scans/lpa-${pop}_predictions/ba_introgressed_095.bed \
    > data/introgression_scans/lpa-${pop}_predictions/ba_introgressed_095.merged.bed

bedtools merge \
    -i data/introgression_scans/lpa-${pop}_predictions/bi_introgressed_095.bed \
    > data/introgression_scans/lpa-${pop}_predictions/bi_introgressed_095.merged.bed
```

##############################################################################
-- 
### Plot introgression probabilities on chromosomes

```

```

### Total size of windows with signals of introgression
```
pop=wel
pop=eel
pop=sel

awk '{sum += $3 - $2} END {print sum}' data/introgression_scans/lpa-${pop}_predictions/ab_introgressed_095.merged.bed
awk '{sum += $3 - $2} END {print sum}' data/introgression_scans/lpa-${pop}_predictions/ba_introgressed_095.merged.bed
awk '{sum += $3 - $2} END {print sum}' data/introgression_scans/lpa-${pop}_predictions/bi_introgressed_095.merged.bed

# wel ab = 166847030 / 2285572469 0.07300010490282118
# wel ba = 210690912 / 2285572469 0.09218299347654595
# wel bi = 12551816 / 2285572469 0.00549176023523409

# eel ab = 802575019 / 2285572469 0.3511483577465161
# eel ba = 180822320 / 2285572469 0.07911467365509293
# eel bi = 13741216 / 2285572469 0.006012155023031125

# sel ab = 332352367 / 2285572469 0.1454131826961554
# sel ba = 148605256 / 2285572469 0.06501883358133852
# sel bi = 5604069 / 2285572469 0.0024519323171808822
```
### Get genes and enrichment:

```
pop=wel
genes_bed=data/genes.bed
bedtools=/Users/enrico/Documents/software/bedtools2/bin/bedtools

for mig in ab ba bi; do
    intro_bed=data/introgression_scans/lpa-${pop}_predictions/${mig}_introgressed_095.merged.bed
    $bedtools intersect \
        -a ${genes_bed} \
        -b ${intro_bed} |
        cut -f4 | cut -d';' -f1 | cut -d'=' -f2 | uniq > data/introgression_scans/lpa-${pop}_genes/${mig}.geneids.txt
done
```
Get table with two columns (gene - GO terms) from gff3:
```
python src/introgression_scans/write_genego_table.py \
    --igff3 'data/LYRU2_2A.FA.gff3' \
    --otable 'data/LYRU2_2A.FA.genego_table.tsv'
```
Enrichment by Lorena:
```
```

### Overlap with genes: random vs observed

For each introgressed window take a random window of the same size from the genome and create a new bed:
```
# create genome file
data/mLynRuf2.2.revcomp.scaffolds.fa.fai | cut -f1,2 | grep "Chr" | grep -v "ChrY" | grep -v "ChrX" \
    > data/mLynRuf2.2.revcomp.scaffolds.genome 

genome=data/mLynRuf2.2.revcomp.scaffolds.genome
bedtools=/Users/enrico/Documents/software/bedtools2/bin/bedtools
pop=wel

for mig in ab ba bi; do
    intro_bed=data/introgression_scans/lpa-${pop}_predictions/${mig}_introgressed_095.merged.bed
    for n in {0..99}; do
        echo "${mig}: ${n}"
        out_bed=data/introgression_scans/lpa-${pop}_randomwins/${mig}/random_windows.${n}.bed
        rm ${out_bed}
        while IFS= read -r line; do
            wsize=$(awk '{sum += $3 - $2} END {print sum}' <(echo "$line"))
            $bedtools random -g ${genome} -n 1 -l $wsize >> ${out_bed}
        done < "$intro_bed"
    done
done
```

total genome = 2285572469
total genome overlap with genes = 937872656
total genes = 937872656 (makes sense)
proportion of genome in genes = 937872656 / 2285572469 = 0.4103447468

Overlap with introgressed windows:

```
bedtools=/Users/enrico/Documents/software/bedtools2/bin/bedtools
pop=wel
genes_bed=data/genes.bed

for mig in ab ba bi; do
    echo $mig > $mig.tmp
    intro_bed=data/introgression_scans/lpa-${pop}_predictions/${mig}_introgressed_095.merged.bed
    $bedtools intersect -a ${intro_bed} -b ${genes_bed} -wo | cut -f8 | awk '{sum += $1} END {print sum}' >> $mig.tmp
    for n in {0..99}; do
        intro_bed=data/introgression_scans/lpa-${pop}_randomwins/${mig}/random_windows.${n}.bed
        $bedtools intersect -a <(cut -f1,2,3 ${intro_bed}) -b ${genes_bed} -wo | cut -f8 | awk '{sum += $1} END {print sum}' >> $mig.tmp
    done
done

paste ab.tmp ba.tmp bi.tmp > data/introgression_scans/lpa-${pop}_randomwins/gene_overlap.table.tsv
rm *.tmp
```

### Length of intro segments distribution

```
Rscript src/introgression_scans/draw_introsegment_lengths.R \
    data/introgression_scans/lpa-wel_predictions \
    plots/introgression_scans/lpa-wel_predictions
```

### pi vs introsegs

Calculate pi along the genome in each population:
```
# vcfdir
vcfdir=/GRUPOS/grupolince/mLynRuf2.2_ref_vcfs
# refdir
refdir=/GRUPOS/grupolince/reference_genomes/lynx_rufus_mLynRuf2.2
# chrs
chrs=($(cat ${refdir}/autosomic_scaffolds_list.txt))

# species1 - eurasian lynx
for pop in wel sel; do
    # vcf
    ivcf=${vcfdir}/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-${pop}.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf
    pop1=data/${pop}.list
    
    # divide vcf per population
    vcftools --vcf $ivcf --keep $pop1 --maf 0.0001 \
        --recode --recode-INFO-all --out data/introgression_scans/${pop}
    
    # calculate pi per population
    vcftools --vcf data/introgression_scans/${pop}.recode.vcf \
        --window-pi 100000 --out data/introgression_scans/${pop}
done

# species2 - iberian lynx
pop2=data/lpa.list

# divide vcf per population
vcftools --vcf $ivcf --keep $pop2 --maf 0.0001 \
    --recode --recode-INFO-all --out data/introgression_scans/lpa
    
# calculate pi per population
vcftools --vcf data/introgression_scans/lpa.recode.vcf \
    --window-pi 100000 --out data/introgression_scans/lpa
```


```
```
