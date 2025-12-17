# Analyses associated with first round of revision

## More evaluation!

Model performance is evaluated in the face of possible genotyping errors (emulating low-depth data) and independent changes in recombination and mutation rates (emulating genome heterogeneity).

### add variable recombination and mutation rates

We aim to emulate genome heterogeneity in mutation and recombination rates by sampling random values around the genome-wide averages of these rates. The [simulate_data_add_murho_change.py](src/simulate_training_data/simulate_data_add_murho_change.py) script ...

### add simulation of low depth 

We aim to emulate the effects of low sequencing data depth by modifying the output of simulations in the step previous to their formatting. The [add_gterror_to_sims.py](src/simulate_training_data/add_gterror_to_sims.py) script ...

```
python src/simulate_training_data/add_gterror_to_sims.py \
    --idir testing_sims_withM/lpa-wel_ab_sims/6_2_5000_fil/sims \
    --odir testing_sims_withM/lpa-wel_ab_sims/6_2_5000_gterror/sims \
    --pop_sizes 40,44 \
    --nbad 12,12 \
    --avg_depth 4
```

### Evaluation data generation:

To evaluate the model on additional simulations, I simulate data under each of the demographic models and each introgression scenario 1000 times (using sampled values for mutation and recombination rates), filter the simulations to 128 SNPs, add genotype errors induced by low-depth in some individuals and generate NPZ files to run the discriminator model on:
```
mkdir data/introgression_scans/revision1/
conda activate lynxtrogression_v2

# for pop_pair in lpa-wel lpa-sel; do
for pop_pair in lpa-wel; do    
    if [ ${pop_pair} == 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
        pop_sizes="40,44"
        n_bad="12,12"
    elif [ ${pop_pair} == 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
        pop_sizes="38,44"
    elif [ ${pop_pair} == 'lpa-sel' ]; then
        models=(12_6 18_7 18_10)
        pop_sizes="24,44"
        n_bad="12,1"
    fi
    for model in ${models[@]}; do
        for migration in ab ba abba baab none; do
            
            # SIM DATA WITH SAMPLED MU AND RHO
            mkdir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}
            python src/introgression_scans/simulate_data_withM_samplemurec.py \
                --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
                --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
                --path_to_msmodified src/introNets/msmodified/ms \
                --migration ${migration} \
                --nreps 2500 \
                --odir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}
            
            # FILTER SIMS
            mkdir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/filtered
            if [ ${migration} == 'ab' ] || [ ${migration} == 'ba' ] || [ ${migration} == 'none' ]; then
                mig=${migration}
                nsims=1000
            elif [ ${migration} == 'abba' ] || [ ${migration} == 'baab' ]; then
                mig="bi"
                nsims=500
            fi
            python src/introgression_scans/filter_sims.py \
                --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration} \
                --odir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/filtered \
                --n_sites 128 --n_sims ${nsims} --migration ${mig} --pop_sizes ${pop_sizes}
            
            # ADD NOISE FROM LOW-DEPTH
            mkdir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/gterror
            python src/introgression_scans/add_gterror_to_sims.py \
                --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/filtered \
                --odir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/gterror \
                --pop_sizes ${pop_sizes} \
                --nbad ${n_bad} \
                --avg_depth 4

            # CREATE NPZ of FILTERED
            python src/introgression_scans/get_npz_from_sims.py \
                --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/filtered \
                --ofile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}.filtered.npz \
                --pop_sizes ${pop_sizes}
            
            # CREATE NPZ of NOISY
            python src/introgression_scans/get_npz_from_sims.py \
                --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/gterror \
                --ofile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}.gterror.npz \
                --pop_sizes ${pop_sizes}
        done
    done
done
```

### Get predictions:

Get predictions from these simulations:
```
conda activate ~/introNets/intronets

# for pop_pair in lpa-wel lpa-sel; do
for pop_pair in lpa-wel; do    
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
        for migration in ab ba abba baab none; do
            for type in filtered gterror; do
                taskset -c 1 python src/introNets/src/models/apply_disc_to_npz.py \
                    --ifile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}.${type}.npz \
                    --ofile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}.${type}.predictions.csv \
                    --weights data/introgression_scans/${pop_pair}_discriminator_withM/test.weights \
                    --pop_sizes ${pop_sizes} \
                    --shape 2,44,128 \
                    --step_size 128 \
                    --in_channels 2 \
                    --n_classes 4
            done
        done
    done
done
```
