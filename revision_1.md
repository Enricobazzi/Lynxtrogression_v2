# Analyses associated with first round of revision

## More evaluation: more realistic evaluation set

Model performance is evaluated in the face of possible genotyping errors (emulating low-depth data) and independent changes in recombination and mutation rates (emulating genome heterogeneity).

In this section I aim to have a more realistic evaluation set. In next sections (see below) I will test for the direct of depth and of recombination and mutation rates.

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
        n_bad="1,12"
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

### plot!

To plot them I run [plot_eval_intro_binary.py](src/introgression_scans/plot_eval_intro_binary.py) and [plot_eval_intro_direction.py](src/introgression_scans/plot_eval_intro_direction.py). I first move all newly generated predictions to separate folders based on evaluation (`data/introgression_scans/revision1/filtered` = just altered rec and mu, `data/introgression_scans/revision1/gterror` = also added low-depth induced genotyping errors):

```
for pop_pair in lpa-wel; do

    if [ ${pop_pair} = 'lpa-wel' ]; then
        models=(12_9 6_2 20_7)
    elif [ ${pop_pair} = 'lpa-eel' ]; then
        models=(34_7 38_4 30_1)
    elif [ ${pop_pair} = 'lpa-sel' ]; then
        models=(12_6 18_7 18_10)
    fi

    for p in '0.9'; do
        pt=$(echo ${p} | tr '.' '_')

        for eval in filtered gterror; do
            echo "plot binary introgression of ${eval} with p_thresh ${p}"
            
            python src/introgression_scans/plot_eval_intro_binary.py \
                --idir data/introgression_scans/revision1/${eval}/ \
                --pop_pair ${pop_pair} \
                --models $(for model in ${models[@]}; do echo ${model}; done | tr '\n' ',') \
                --pthresh ${p} \
                --oplot plots/introgression_scans/revision1/${pop_pair}.intro_binary.${eval}.${pt}.cm.pdf
            
            echo "plot introgression direction of ${eval} with p_tresh ${p}"
            python src/introgression_scans/plot_eval_intro_direction.py \
                --idir data/introgression_scans/revision1/${eval}/ \
                --pop_pair ${pop_pair} \
                --models $(for model in ${models[@]}; do echo ${model}; done | tr '\n' ',') \
                --pthresh ${p} \
                --oplot plots/introgression_scans/revision1/${pop_pair}.intro_direction.${eval}.${pt}.cm.pdf
        done
    done
done
```
To get precision - recall graphs I run [plot_eval_intro_binary.py](src/introgression_scans/plot_eval_intro_binary.py) 

## More evaluation: check 2X to 6X depth in all samples

To more directly test for the effect of low sequencing depth on model performance I take the same simulations as above but add gt errors assuming all samples have 2X to 6X.

### add GT errors at different depths

I use the same simulations as above:
```
conda activate lynxtrogression_v2

for dp in {2..6}; do
    for pop_pair in lpa-wel; do    
        if [ ${pop_pair} == 'lpa-wel' ]; then
            models=(12_9 6_2 20_7)
            pop_sizes="40,44"
            n_bad="20,22"
        elif [ ${pop_pair} == 'lpa-eel' ]; then
            models=(34_7 38_4 30_1)
            pop_sizes="38,44"
        elif [ ${pop_pair} == 'lpa-sel' ]; then
            models=(12_6 18_7 18_10)
            pop_sizes="24,44"
            n_bad="12,22"
        fi
        for model in ${models[@]}; do
            for migration in ab ba abba baab none; do
                echo "${pop_pair} ${model} ${migration} ${dp}X"
                # ADD NOISE FROM LOW-DEPTH
                mkdir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/gterror_${dp}X
                python src/introgression_scans/add_gterror_to_sims.py \
                    --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/filtered \
                    --odir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/gterror_${dp}X \
                    --pop_sizes ${pop_sizes} \
                    --nbad ${n_bad} \
                    --avg_depth ${dp}
                
                # CREATE NPZ of NOISY
                python src/introgression_scans/get_npz_from_sims.py \
                    --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}/gterror_${dp}X \
                    --ofile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}.gterror_${dp}X.npz \
                    --pop_sizes ${pop_sizes}
            done
        done
    done
done
```

### Get predictions at different depths

```
conda activate ~/introNets/intronets

for dp in {2..6}; do
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
                taskset -c 1 python src/introNets/src/models/apply_disc_to_npz.py \
                    --ifile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}.gterror_${dp}X.npz \
                    --ofile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}.gterror_${dp}X.predictions.csv \
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

...

## More evaluation: check the effects of extreme mu and rec rates

To more directly test for the effect of extremely high or low mutation and recombination rates I repeat simulations with the following combinations:
- mu / 10 & rec / 10
- mu / 10 & rec * 10
- mu * 10 & rec / 10
- mu * 10 & rec * 10

### Simulate data with extreme mu and rec

```
conda activate lynxtrogression_v2

for mu in 6e-8 6e-10; do
    for rec in 1.9e-7 1.9e-9; do
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
                n_bad="1,12"
            fi
            for model in ${models[@]}; do
                for migration in ab ba abba baab none; do
                    
                    # SIM DATA WITH SAMPLED MU AND RHO
                    mkdir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}_mu${mu}_rec${rec}
                    python src/introgression_scans/simulate_data_withM_samplemurec.py \
                        --demes_yaml data/demographic_inference/${pop_pair}_best_yamls/${pop_pair}_${model}_final_best_model.yaml \
                        --confint data/demographic_inference/${pop_pair}_CI/${pop_pair}.${model}.CI.csv \
                        --path_to_msmodified src/introNets/msmodified/ms \
                        --migration ${migration} \
                        --nreps 2500 \
                        --mu ${mu} \
                        --rec ${rec} \
                        --odir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}_mu${mu}_rec${rec}
                    
                    # FILTER SIMS
                    mkdir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}_mu${mu}_rec${rec}/filtered
                    if [ ${migration} == 'ab' ] || [ ${migration} == 'ba' ] || [ ${migration} == 'none' ]; then
                        mig=${migration}
                        nsims=1000
                    elif [ ${migration} == 'abba' ] || [ ${migration} == 'baab' ]; then
                        mig="bi"
                        nsims=500
                    fi
                    python src/introgression_scans/filter_sims.py \
                        --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}_mu${mu}_rec${rec} \
                        --odir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}_mu${mu}_rec${rec}/filtered \
                        --n_sites 128 --n_sims ${nsims} --migration ${mig} --pop_sizes ${pop_sizes}
                    
                    # CREATE NPZ of FILTERED
                    python src/introgression_scans/get_npz_from_sims.py \
                        --idir data/introgression_scans/revision1/${pop_pair}_${model}_${migration}_mu${mu}_rec${rec}/filtered \
                        --ofile data/introgression_scans/revision1/${pop_pair}_${model}_${migration}_mu${mu}_rec${rec}.filtered.npz \
                        --pop_sizes ${pop_sizes}
                    
                done
            done
        done
    done
done
```
