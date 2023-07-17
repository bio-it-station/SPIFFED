echo ">> Running PRINCE_DATA!!"
for sample in "condition1" "condition2"
do
    echo "python ./main.py -s 11101001 /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/gold_standard.tsv /ccb/salz3/kh.chao/SPIFFED/output/EPIC_res/PRINCE_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M RF -n 6 -m EXP -f STRING --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 5"

    python ./main.py -s 11101001 /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/gold_standard.tsv       ${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M RF -n 6 -m EXP -f STRING --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 5
done
