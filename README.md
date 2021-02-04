# EPPC: the tool for improving EPIC results

Here is the main and only one command that you need to run EPPC.
```
python ./main.py -s `feature_selection` `input directory` -c `path to gold standard file` `output directory` -o `output filename prefix` -M `training method` -n `number of cores` -m EXP -f STRING --LEARNING_SELECTION `learning method selection` --K_D_TRAIN `k-fold or direct training` --FOLD_NUM `number of folds` --TRAIN_TEST_RATIO `testing dataset ratio` --POS_NEG_RATIO `negative PPIs ratio` --NUM_EP `number of elution profiles` --NUM_FRC `number of fractions`
```


## Supervised learning: EPIC scores + random Forest (EPIC scores + RF)

## Supervised Learning: raw inputs + convolutional neural networks (Raw + CNN)

## Semi-supervised learning: raw inputs + convolutional neural networks (Raw + CNN)

## Semi-supervised learning: raw inputs + label spreading (Raw + LS)
