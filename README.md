# EPIC Improvement

## EPIC + RF

### TEST
```
python ./main.py -s 11101001 /home/kuan-hao/EPPC/input/real_data/elution_profiles_test/ -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/TEST/ -o TEST -M RF -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --N_ESTIMATORS 100 --MAX_DEPTH 10 --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 1 --NUM_EP 2 --NUM_FRC 27
```


### SELF DATA
```
python ./main.py -s 11101001 /home/kuan-hao/EPPC/input/real_data/elution_profiles -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/TEST/ -o Intensity_H -M RF -n 10 -m EXP -f STRING --N_ESTIMATORS 100 --MAX_DEPTH 10 --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO
```

### EPIC DATA

```
python ./main.py -s 11101001 /home/kuan-hao/EPPC/input/EPIC_test -c /home/kuan-hao/EPPC/input/Worm_reference_complexes.txt /home/kuan-hao/EPPC/output/TEST/ -o Beads_A -M RF -n 10 -m EXP -f STRING --N_ESTIMATORS 100 --MAX_DEPTH 10 --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 1
```


## RAW + CNN

### TEST (CNN+ssl)
```
python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/real_data/elution_profiles_test/ -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/TEST2/ -o TEST2 -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION ssl --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.7 --POS_NEG_RATIO 1 --NUM_EP 2 --NUM_FRC 27
```

### TEST (CNN+ssl)
```
python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/real_data/elution_profiles_test/ -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/TEST/ -o TEST -M CNN -n 6 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 1 --NUM_EP 2 --NUM_FRC 27 --EVAL_RATIO 0.3
```

### SELF DATA
```
python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/real_data/elution_profiles/ -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/TEST/ -o TEST -M LS -n 6 -m EXP -f STRING --LEARNING_SELECTION ssl --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.7 --POS_NEG_RATIO 1 --NUM_EP 2 --NUM_FRC 27 --EVAL_RATIO 0.3
```

### EPIC DATA
```
python ./main.py -s 000000001 /Users/chaokuan-hao/Documents/BIO_IT_Station/input/EPIC_test -c /Users/chaokuan-hao/Documents/BIO_IT_Station/input/Worm_reference_complexes.txt /Users/chaokuan-hao/Documents/BIO_IT_Station/output/EPIC_DATA_RESULT/RAW_CNN_DIRECT_TRAIN/1_5_w/30_70/RAW_1_CNN_ratio_1_1_BEADS_A_Conv2D_Dense_60_30_Drp_50_20_Drp_50_10_5_1 -o RAW_1_CNN_ratio_1_1_BEADS_A_Conv2D_Dense_60_30_Drp_50_20_Drp_50_10_5_1 -M CNN -n 6 -m EXP -f STRING
```

## RAW + CNN (ENSEMBLE)

### TEST
```
python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/OUR_DATA/intensity_HML_ensemble/ -c /home/kuan-hao/EPPC/input/OUR_DATA/gold_standard.tsv /home/kuan-hao/EPPC/output/SELF_DATA/intensity_HML_ensemble__negative_ratio_5/ -o out -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.7 --POS_NEG_RATIO 5 --NUM_EP 2 --NUM_FRC 27 --CNN_ENSEMBLE 1
```
```
python ./src/main.py -s 000000001 /home/kuan-hao/EPPC/input/EPIC_DATA/beads_ensemble/ -c /home/kuan-hao/EPPC/input/EPIC_DATA/Worm_reference_complexes.txt /home/kuan-hao/EPPC/output/EPIC_DATA/beads_ensemble__negative_ratio_5/ -o out -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.7 --POS_NEG_RATIO 5 --NUM_EP 2 --NUM_FRC 60 --CNN_ENSEMBLE 1

```

### TEST2
```
python ./main.py -s 000000001 /Users/chaokuan-hao/Documents/BIO_IT_Station/input/EPIC_Multi -c /Users/chaokuan-hao/Documents/BIO_IT_Station/input/Worm_reference_complexes.txt /Users/chaokuan-hao/Documents/BIO_IT_Station/output/EPIC_DATA_RESULT/RAW_CNN_ENSEMBLE/1_1/20_80/ -o BEADS -M CNN -n 6 -m EXP -f STRING
```

### EPIC Data
```
python ./main.py -s 000000001 /Users/chaokuan-hao/Documents/BIO_IT_Station/input/EPIC_Multi/ -c /Users/chaokuan-hao/Documents/BIO_IT_Station/input/Worm_reference_complexes.txt /Users/chaokuan-hao/Documents/BIO_IT_Station/output/EPIC_DATA_RESULT/RAW_CNN_ENSEMBLE/1_1/20_80/ -o BEADS -M CNN -n 6 -m EXP -f STRING
```

## RAW + LSTM
```
python ./main.py -s 000000001 /Users/chaokuan-hao/Documents/BIO_IT_Station/input/real_data/elution_profiles -c /Users/chaokuan-hao/Documents/BIO_IT_Station/input/real_data/gold_standard.tsv /Users/chaokuan-hao/Documents/BIO_IT_Station/output/SELF_DATA_RESULT/RAW_1_LSTM_intensity_H_LSTM_128_Dr_20_LSTM_128_Dr_20_Dense_1_POS_CUTOFF_30 -o RAW_1_LSTM_intensity_H_LSTM_128_Dr_20_LSTM_128_Dr_20_Dense_1_POS_CUTOFF_30 -M LSTM -n 6 -m EXP -f STRING
```





```
python ./lcl_cutoff.py -s 000000001 /Users/chaokuan-hao/Documents/BIO_IT_Station/input/real_data/elution_profiles_ensemble/intensity_L -c /Users/chaokuan-hao/Documents/BIO_IT_Station/input/real_data/gold_standard.tsv  /Users/chaokuan-hao/Documents/BIO_IT_Station/output/SELF_DATA_RESULT/RAW_CNN_ENSEMBLE/1_1/20_80/intensity_L -o intensity -M CNN -n 6 -m EXP -f STRING
```
