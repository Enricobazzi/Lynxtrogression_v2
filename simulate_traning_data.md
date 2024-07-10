## Simulating Training Data

### Download Intronets

I need [introNets](https://github.com/SchriderLab/introNets) version of ms (msmodified) to run the simulations. I download the folder from github and compile msmodified as suggested by Dylan.
```
cd src
git clone https://github.com/SchriderLab/introNets.git
cd introNets/msmodified
```

Script that takes a model from gadma yaml and confidence param table from gadma-get_confidence_intervals and spits ms command custom migration (no, 1->2, 2->1, 1->2 and 2->1)

I manually modify the output demes yamls from gadma2 that have 4 populations (wel, lpa, wel-lpa, ancestral) to only include 2 populations (wel, lpa) so that the msmodified command is the same as in introunets.

```
pop_pair=lpa-wel
for model in 12_9 20_7 6_2; do
    for migration in lpa-lly lly-lpa none bi; do

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
            --nreps 10000 \
            --odir data/simulate_training_data/${pop_pair}_${migration}_sims/${model}/
    
    done
done

for model in 12_9 20_7 6_2; do
    for migration in bi; do

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
            --nreps 10000 \
            --odir data/simulate_training_data/${pop_pair}_${migration}_sims/${model}/
    
    done
done

```


### format training data for training

```
conda activate ~/introNets/intronets

pop_pair=lpa-wel
migration=bi

mpirun -n 4 python src/introNets/src/data/format.py \
    --verbose \
    --idir data/simulate_training_data/${pop_pair}_${migration}_sims/ \
    --ofile test_${pop_pair}_${migration}.hdf5 \
    --pop_sizes 40,44 --out_shape 2,84,128 --pop -1 \
    --sorting seriate_match --metric cosine
```