screen -S 2_ChrA1
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrA1
taskset -c 1 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrC1
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrC1
taskset -c 2 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrB1
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrB1
taskset -c 3 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrA2_rc
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrA2_rc
taskset -c 4 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrC2
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrC2
taskset -c 5 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrB2_rc
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrB2_rc
taskset -c 18 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrB3
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrB3
taskset -c 6 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrB4_rc
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrB4_rc
taskset -c 7 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrA3_rc
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrA3_rc
taskset -c 8 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrD1
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrD1
taskset -c 9 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrD4
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrD4
taskset -c 10 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrD3
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrD3
taskset -c 11 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrD2
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrD2
taskset -c 12 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrF2
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrF2
taskset -c 13 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrF1_rc
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrF1_rc
taskset -c 14 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrE2_rc
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrE2_rc
taskset -c 15 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrE1
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrE1
taskset -c 16 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4

screen -S 2_ChrE3_rc
conda activate ~/introNets/intronets
pop="wel"
pop_sizes="40,44"
chr=mLynRuf2.2_ChrE3_rc
taskset -c 17 python src/introNets/src/models/apply_disc_to_npz.py \
    --ifile data/introgression_scans/npz_files/lpa-${pop}.${chr}.npz \
    --ofile data/introgression_scans/lpa-${pop}_predictions_withM/${chr}.predictions.csv \
    --weights data/introgression_scans/lpa-${pop}_discriminator_withM/test.weights \
    --pop_sizes ${pop_sizes} \
    --shape 2,44,128 \
    --step_size 64 \
    --in_channels 2 \
    --n_classes 4
