# echo ">> Running EPIC_DATA!!"
# for sample in "beadsA" "beadsALF" "beadsB" "beadsBNF"
# do
#     echo "python ./main.py -s 000000001 /ccb/salz3/kh.chao/SPIFFED/input/EPIC_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/EPIC_DATA/Worm_reference_complexes.txt /ccb/salz3/kh.chao/SPIFFED/output/EPIC_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 50 --CNN_ENSEMBLE 0"

#     python ./main.py -s 000000001 /ccb/salz3/kh.chao/SPIFFED/input/EPIC_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/EPIC_DATA/Worm_reference_complexes.txt /ccb/salz3/kh.chao/SPIFFED/output/EPIC_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 50 --CNN_ENSEMBLE 0
# done

echo ">> Running OUR_DATA!!"
# for sample in "intensity_high" "intensity_medium" "intensity_low"
for sample in "intensity_medium" "intensity_low"
do
    echo "python ./main.py -s 000000001 /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/gold_standard.tsv /ccb/salz3/kh.chao/SPIFFED/output/OUR_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 50 --CNN_ENSEMBLE 0"

    python ./main.py -s 000000001 /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/gold_standard.tsv /ccb/salz3/kh.chao/SPIFFED/output/OUR_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 50 --CNN_ENSEMBLE 0
done