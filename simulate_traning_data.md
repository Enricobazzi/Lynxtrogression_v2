## Simulating Training Data

### Download Intronets

I need [introNets](https://github.com/SchriderLab/introNets) version of ms (msmodified) to run the simulations. I download the folder from github and compile msmodified as suggested by Dylan.
```
cd src
git clone https://github.com/SchriderLab/introNets.git
cd introNets/msmodified
```

### Run simulations

Script that takes a model from gadma yaml and confidence param table from gadma-get_confidence_intervals and spits ms command custom migration (no, 1->2, 2->1, 1->2 and 2->1)

```
pop_pair=lpa-wel
for model in 12_9 20_7 6_2; do
    for migration in ab ba none bi; do

        if [ -d data/simulate_training_data/${pop_pair}_${migration}_sims/${model} ]; then
            rm -rf data/simulate_training_data/${pop_pair}_${migration}_sims/${model}/*
        else
            mkdir -p data/simulate_training_data/${pop_pair}_${migration}_sims/${model}/
        fi
        
        echo "simulating ${pop_pair}_${migration} of ${model}"

        python src/simulate_training_data/simulate_data.py \
            --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
            --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
            --path_to_msmodified src/introNets/msmodified/ms \
            --migration ${migration} \
            --nreps 15000 \
            --odir data/simulate_training_data/${pop_pair}_${migration}_sims/${model}/
    
    done
done

for migration in ab ba none bi; do
    for model in 12_9 20_7 6_2; do
        echo "${pop_pair}_${migration}_${model}"
        zgrep "segsites" data/simulate_training_data/${pop_pair}_${migration}_sims/${model}/mig.msOut.gz |
            cut -d' ' -f2 | awk '$1 > 128' | wc -l
    done
done | grep -v "lpa" | awk '{sum += $1} END {print sum}'
```

### Format training data for training

I filter out simulations with < 128 segsites (or format.py will throw errors and corrupt the hdf5 file):

```
pop_pair=lpa-wel
for model in 12_9 20_7 6_2; do
    for migration in ab ba none bi; do

        if [ -d data/simulate_training_data/${pop_pair}_${migration}_filtered_sims/${model} ]; then
            rm -rf data/simulate_training_data/${pop_pair}_${migration}_filtered_sims/${model}/*
        else
            mkdir -p data/simulate_training_data/${pop_pair}_${migration}_filtered_sims/${model}/
        fi
        
        echo "filtering ${pop_pair}_${migration} of ${model}"

        python src/simulate_training_data/filter_sims.py \
            --idir data/simulate_training_data/${pop_pair}_${migration}_sims/${model}/ \
            --odir data/simulate_training_data/${pop_pair}_${migration}_filtered_sims/${model}/ \
            --n_sites 128 --n_sims 10000
    done
done
```

To format filtered simulation:
```
conda activate ~/introNets/intronets
pop_pair=lpa-wel

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_ab_filtered_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_ab.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 1 |& tee format_ab.log

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_ba_filtered_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_ba.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 0 |& tee format_ba.log

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_bi_filtered_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_bi.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop -1 |& tee format_bi.log

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_none_filtered_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_none.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --include_zeros |& tee format_none.log
```

Check hdf5:
```
conda activate ~/introNets/intronets
pop_pair=lpa-wel
for migration in ab ba bi none; do
    python src/introNets/src/data/get_h5_stats.py \
        --ifile data/simulate_training_data/${pop_pair}_${migration}.hdf5
done
```

### Train disc model

```
script train_disc.log

conda activate ~/introNets/intronets
pop_pair=lpa-wel
taskset -c 1,2,3,4,5,6,7,8,9,10 \
    python src/introNets/src/models/train_discriminator.py \
        --idir data/simulate_training_data/ --odir train_disc --n_classes 4
```

##########

### minitest

```
pop_pair=lpa-wel
conda activate ~/introNets/intronets

for model in 12_9 20_7 6_2; do
    for migration in ab ba none bi; do

        odir=minitest/${migration}/${model}

        if [ -d ${odir} ]; then
            rm -rf ${odir}/*
        else
            mkdir -p ${odir}/
        fi
        
        echo "simulating ${pop_pair}_${migration} of ${model}"

        python src/simulate_training_data/simulate_data.py \
            --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
            --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
            --path_to_msmodified src/introNets/msmodified/ms \
            --migration ${migration} \
            --nreps 100 \
            --odir ${odir}/
    
    done
done

for migration in ab ba none bi; do
    for model in 12_9 20_7 6_2; do
        echo "${migration} ${model}"
        zgrep "segsites" minitest/${migration}/${model}/mig.msOut.gz |
            cut -d' ' -f2 | awk '$1 > 128' | wc -l
    done
done | grep -v "lpa" | awk '{sum += $1} END {print sum}'

for model in 12_9 20_7 6_2; do
    for migration in ab ba none bi; do

        if [ -d minitest/${migration}_fil/${model} ]; then
            rm -rf minitest/${migration}_fil/${model}/*
        else
            mkdir -p minitest/${migration}_fil/${model}
        fi
        
        echo "filtering ${pop_pair}_${migration} of ${model}"

        python src/simulate_training_data/filter_sims.py \
            --idir minitest/${migration}/${model}/ \
            --odir minitest/${migration}_fil/${model}/ \
            --n_sites 128 --n_sims 60
    done
done

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir minitest/ab_fil/ \
    --ofile minitest/${pop_pair}_ab.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 1 --chunk_size 2

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir minitest/ba_fil/ \
    --ofile minitest/${pop_pair}_ba.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 0 --chunk_size 2

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir minitest/bi_fil/ \
    --ofile minitest/${pop_pair}_bi.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop -1 --chunk_size 2

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir minitest/none_fil/ \
    --ofile minitest/${pop_pair}_none.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --include_zeros --chunk_size 2
```

### real data

VCF to NPZ
```
vcfFileName=data/phased_data/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-wel.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf
species1ListFileName=data/phased_data/wel.list
species2ListFileName=data/phased_data/lpa.list
maskedRefFileName=data/mLynRuf2.2.revcomp.scaffolds.fa
arm=mLynRuf2.2_ChrA1
npzFileName=test.npz

python phasedVcfToNpz.py $vcfFileName $species1ListFileName $species2ListFileName $maskedRefFileName $arm $npzFileName
```

NPZ to HDF5
```
conda activate ~/introNets/intronets

mpirun -n 8 python src/introNets/src/data/format_npz.py \
    --verbose \
    --ifile test.npz \
    --ofile test.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 \
    --keys sechMatrix,simMatrix
```

```
cat intro.txt | cut -d' ' -f2,4,6,7,8,9 | tr ':' ' ' | tr ' ' '\t' | sort -nk 1| tr '\t' ' ' | grep -n "," | tr -d "\'" | tr -d "[" | tr -d "]" | tr -d ',' | tr ' ' '\t' | tr ':' '\t' > intro_table.tsv
```