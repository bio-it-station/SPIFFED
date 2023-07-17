echo ">> Running OUR_DATA!!"
for sample in "intensity_high" "intensity_medium" "intensity_low"
do
    echo "python ./main.py -s 11101001 /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/gold_standard.tsv /ccb/salz3/kh.chao/SPIFFED/output/EPIC_res/OUR_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_50__test_ratio_30 -o TEST -M RF -n 6 -m EXP -f STRING --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 50"

    python ./main.py -s 11101001 /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/OUR_DATA/gold_standard.tsv /ccb/salz3/kh.chao/SPIFFED/output/EPIC_res/OUR_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_50__test_ratio_30 -o TEST -M RF -n 6 -m EXP -f STRING --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 50
done
