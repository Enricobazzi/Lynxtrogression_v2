## Simulating Training Data

### Download Intronets

I need [introNets](https://github.com/SchriderLab/introNets) version of ms (msmodified) to run the simulations. I download the folder from github and compile msmodified as suggested by Dylan.
```
cd src
git clone https://github.com/SchriderLab/introNets.git
cd introNets/msmodified
```

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
        zgrep "segsites" ${pop_pair}_${migration}_sims/${model}/mig.msOut.gz |
            cut -d' ' -f2 | awk '$1 > 128' | wc -l
    done
done | grep -v "lpa" | awk '{sum += $1} END {print sum}'
```

### format training data for training

```
script format_ab.log
conda activate ~/introNets/intronets
pop_pair=lpa-wel
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_ab_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_ab.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 1 \
    --sorting seriate_match --metric cosine

script format_ba.log
conda activate ~/introNets/intronets
pop_pair=lpa-wel
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_ba_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_ba.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 0 \
    --sorting seriate_match --metric cosine

script format_bi.log
conda activate ~/introNets/intronets
pop_pair=lpa-wel
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_bi_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_bi.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop -1 \
    --sorting seriate_match --metric cosine

script format_none.log
conda activate ~/introNets/intronets
pop_pair=lpa-wel
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_none_sims/ \
    --ofile data/simulate_training_data/${pop_pair}_none.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --include_zeros \
    --sorting seriate_match --metric cosine
```

```
conda activate ~/introNets/intronets
pop_pair=lpa-wel
for migration in ab ba bi none; do
    python src/introNets/src/data/get_h5_stats.py \
        --ifile data/simulate_training_data/${pop_pair}_${migration}.hdf5
done
```

```
script train_disc.log

conda activate ~/introNets/intronets
pop_pair=lpa-wel
taskset -c 1,2,3,4,5,6,7,8,9,10 \
    python src/introNets/src/models/train_discriminator.py \
        --idir data/simulate_training_data/ --odir train_disc --n_classes 4 --batch_size 16

6+1
```

```
script train_disc.log

module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eye/eba/intronets

pop_pair=lpa-wel
python src/introNets/src/models/train_discriminator.py \
    --idir data/simulate_training_data/ --odir train_disc --n_classes 4 --batch_size 16
```

##########

script test.log
conda activate ~/introNets/intronets
pop_pair=lpa-wel
mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_ab_sims/ \
    --ofile test.${pop_pair}_ab.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 1 \
    --sorting seriate_match --metric cosine --chunk_size 16

script test2.log
conda activate ~/introNets/intronets
pop_pair=lpa-wel
mpirun -n 12 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_ab_sims/ \
    --ofile test2.${pop_pair}_ab.hdf5 \
    --pop_sizes 40,44 --out_shape 2,44,128 --pop 1 \
    --sorting seriate_match --metric cosine
