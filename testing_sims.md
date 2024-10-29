## Testing sims

I test how simulated data and trained model look with different timings of migration in the different population pairs:

```
mkdir testing_sims
conda activate lynxtrogression_v2
```

### Lowering migration time

I modify [simulate_data.py](src/introgression_scans/simulate_data.py) so that the migration time is a up to 5k, 10k or 20k generations before present.

```
migration=ab

for pop_pair in lpa-wel lpa-eel; do
    mkdir testing_sims/${pop_pair}_${migration}_sims/
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    fi
    for model in ${models[@]}; do
        for mtime in 5000 10000 20000; do
            mkdir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}
            echo "testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}/"
            python src/introgression_scans/simulate_data.py \
                --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
                --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
                --path_to_msmodified src/introNets/msmodified/ms \
                --migration ${migration} \
                --nreps 25 \
                --odir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime} \
                --mtime ${mtime}
        done
    done
done

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    fi
    for model in ${models[@]}; do
        for mtime in 5000 10000 20000; do
            echo "testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
            mkdir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/
            mkdir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims
            python src/introgression_scans/filter_sims.py \
                --idir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}/ \
                --odir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims \
                --n_sites 128 --n_sims 15
        done
    done
done
```

Write the hdf5
```
conda activate ~/introNets/intronets
migration=ab

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    fi
    for model in ${models[@]}; do
        for mtime in 5000 10000 20000; do
            echo "testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
            mpirun -n 4 python src/introNets/src/data/format.py \
                --verbose \
                --idir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
                --ofile testing_sims/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
                --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 1
        done
    done
done
```

Write fasta
```
conda activate ~/introNets/intronets
migration=ab
mkdir testing_sims/fastas/

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    fi
    for model in ${models[@]}; do
        for mtime in 5000 10000 20000; do
            echo "testing_sims/fastas/${pop_pair}_${model}_${mtime}_${migration}"
            mkdir testing_sims/fastas/${pop_pair}_${model}_${mtime}_${migration}
            python src/introgression_scans/write_fastas_from_hdf5.py \
                --ifile testing_sims/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
                --odir testing_sims/fastas/${pop_pair}_${model}_${mtime}_${migration} \
                --nseqs 15
        done
    done
done
```

### Checking all intro scenarios

SIM DATA
```
mtime=10000

for pop_pair in lpa-wel lpa-eel lpa-sel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    elif [ ${pop_pair} == 'lpa-sel' ]; then
        models=(12_6 18_7 18_10)
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            mkdir testing_sims/${pop_pair}_${migration}_sims/
            mkdir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}
            echo "testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}/"
            python src/introgression_scans/simulate_data.py \
                --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
                --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
                --path_to_msmodified src/introNets/msmodified/ms \
                --migration ${migration} \
                --nreps 30 \
                --odir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime} \
                --mtime ${mtime}
        done
    done
done
```

FILTER SIMS
```
mtime=5000
for pop_pair in lpa-wel lpa-eel lpa-sel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    elif [ ${pop_pair} == 'lpa-sel' ]; then
        models=(12_6 18_7 18_10)
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            echo "testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
            mkdir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/
            mkdir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims
            python src/introgression_scans/filter_sims.py \
                --idir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}/ \
                --odir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims \
                --n_sites 128 --n_sims 20
        done
    done
done
```

MAKE HDF5
```
conda activate ~/introNets/intronets
mtime=10000

for pop_pair in lpa-wel lpa-eel lpa-sel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    elif [ ${pop_pair} == 'lpa-sel' ]; then
        models=(12_6 18_7 18_10)
        pop_sizes="24,44"
    fi
    for model in ${models[@]}; do
        echo "testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
        migration=ab
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 1
        migration=ba
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 0
        migration=none
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --include_zeros
    done
done
```

GET FASTA
```
conda activate ~/introNets/intronets
mtime=10000
mkdir testing_sims/fastas/

for pop_pair in lpa-wel lpa-eel lpa-sel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    elif [ ${pop_pair} == 'lpa-sel' ]; then
        models=(12_6 18_7 18_10)
        pop_sizes="24,44"
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            echo "testing_sims/fastas/${pop_pair}_${model}_${mtime}_${migration}"
            mkdir testing_sims/fastas/${pop_pair}_${model}_${mtime}_${migration}
            python src/introgression_scans/write_fastas_from_hdf5.py \
                --ifile testing_sims/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
                --odir testing_sims/fastas/${pop_pair}_${model}_${mtime}_${migration} \
                --nseqs 20
        done
    done
done
```

### Sim with mig rate

```
mkdir testing_sims_withM
conda activate lynxtrogression_v2
```

SIM DATA
```
mtime=5000

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            mkdir testing_sims_withM/${pop_pair}_${migration}_sims/
            mkdir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}
            echo "testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}/"
            python src/introgression_scans/simulate_data_withM.py \
                --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
                --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
                --path_to_msmodified src/introNets/msmodified/ms \
                --migration ${migration} \
                --nreps 30 \
                --odir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime} \
                --mtime ${mtime}
        done
    done
done
```

FILTER SIMS
```
mtime=5000
for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            echo "testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
            mkdir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/
            mkdir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims
            python src/introgression_scans/filter_sims.py \
                --idir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}/ \
                --odir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims \
                --n_sites 128 --n_sims 15
        done
    done
done
```

MAKE HDF5
```
conda activate ~/introNets/intronets
mtime=5000

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    fi
    for model in ${models[@]}; do
        echo "testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
        migration=ab
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims_withM/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 1
        migration=ba
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims_withM/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 0
        migration=none
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims_withM/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims_withM/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --include_zeros
    done
done
```

GET FASTA
```
conda activate ~/introNets/intronets
mtime=5000
mkdir testing_sims_withM/fastas/

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            echo "testing_sims_withM/fastas/${pop_pair}_${model}_${mtime}_${migration}"
            mkdir testing_sims_withM/fastas/${pop_pair}_${model}_${mtime}_${migration}
            python src/introgression_scans/write_fastas_from_hdf5.py \
                --ifile testing_sims_withM/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
                --odir testing_sims_withM/fastas/${pop_pair}_${model}_${mtime}_${migration} \
                --nseqs 20
        done
    done
done
```

### Sim with fixed mig proportion

```
mkdir testing_sims_fix
conda activate lynxtrogression_v2
```

SIM DATA
```
mtime=10000

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            mkdir testing_sims_fix/${pop_pair}_${migration}_sims/
            mkdir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}
            echo "testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}/"
            python src/introgression_scans/simulate_data_fix.py \
                --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
                --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
                --path_to_msmodified src/introNets/msmodified/ms \
                --migration ${migration} \
                --nreps 20 \
                --odir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime} \
                --mtime ${mtime}
        done
    done
done
```

FILTER SIMS
```
mtime=10000
for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            echo "testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
            mkdir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/
            mkdir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims
            python src/introgression_scans/filter_sims.py \
                --idir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}/ \
                --odir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/sims \
                --n_sites 128 --n_sims 10
        done
    done
done
```

MAKE HDF5
```
conda activate ~/introNets/intronets
mtime=10000

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    fi
    for model in ${models[@]}; do
        echo "testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/"
        migration=ab
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims_fix/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 1
        migration=ba
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims_fix/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --pop 0
        migration=none
        mpirun -n 4 python src/introNets/src/data/format.py \
            --verbose \
            --idir testing_sims_fix/${pop_pair}_${migration}_sims/${model}_${mtime}_fil/ \
            --ofile testing_sims_fix/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
            --pop_sizes ${pop_sizes} --out_shape 2,44,128 --include_zeros
    done
done
```

GET FASTA
```
conda activate ~/introNets/intronets
mtime=10000
mkdir testing_sims_fix/fastas/

for pop_pair in lpa-wel lpa-eel; do
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    fi
    for model in ${models[@]}; do
        for migration in ab ba none; do
            echo "testing_sims_fix/fastas/${pop_pair}_${model}_${mtime}_${migration}"
            mkdir testing_sims_fix/fastas/${pop_pair}_${model}_${mtime}_${migration}
            python src/introgression_scans/write_fastas_from_hdf5.py \
                --ifile testing_sims_fix/${pop_pair}_${model}_${mtime}_${migration}.hdf5 \
                --odir testing_sims_fix/fastas/${pop_pair}_${model}_${mtime}_${migration} \
                --nseqs 15
        done
    done
done
```
