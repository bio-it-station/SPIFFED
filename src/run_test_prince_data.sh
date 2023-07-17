echo ">> Running OUR_DATA!!"
for sample in "condition1" "condition2"
do
    echo "python ./main.py -s 000000001 /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/gold_standard.tsv /ccb/salz3/kh.chao/SPIFFED/output/PRINCE_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 5 --CNN_ENSEMBLE 0"

    python ./main.py -s 000000001 /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/${sample} -c /ccb/salz3/kh.chao/SPIFFED/input/PRINCE_DATA/gold_standard.tsv /ccb/salz3/kh.chao/SPIFFED/output/PRINCE_DATA/${sample}/TEST/CNN_SL/FOLDS/${sample}__K_D__k__CNN_SL__fold_number_5__n    egative_ratio_5__test_ratio_30 -o TEST -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 5 --CNN_ENSEMBLE 0
done