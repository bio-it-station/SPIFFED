# EPPC: the tool for improving EPIC results

## Main goal
Improve the current best elution profile-based protein complexes inference tool -- EPIC

---

## Installation

---

## Run EPPC
Here is the main and only one command that you need to run:

> python ./main.py -s <b>`feature_selection`</b> <b>`input_directory`</b> -c <b>`gold_standard_file_path`</b> <b>`output_directory`</b> -o <b>`output_filename_prefix`</b> -M <b>`training_method`</b> -n <b>`number_of_cores`</b> -m EXP -f STRING --LEARNING_SELECTION <b>`learning method selection`</b> --K_D_TRAIN <b>`fold_or_direct_training`</b> --FOLD_NUM <b>`number_of_folds`</b> --TRAIN_TEST_RATIO <b>`testing_data_ratio`</b> --POS_NEG_RATIO <b>`negative_PPIs_ratio`</b> --NUM_EP <b>`number_of_elution_profiles`</b> --NUM_FRC <b>`number_of_fractions`</b>

Don't worry if the above command looks a bit scary. In the next section, we detailly explain the definition and possible options for each parameter.

---

## Parameter Definition

1. (`-s` <b>`feature_selection`</b>) or (`--feature_selection` <b>`feature_selection`</b>): Specify correlation scores to be used in EPIC. Eight different correlation socres are implemented in EPIC, in order: <b>Mutual Information</b>, <b>Bayes Correlation</b>, <b>Euclidean Distance</b>, <b>Weighted Cross-Correlation</b>, <b>Jaccard Score</b>, <b>PCCN</b>, <b>Pearson Correlation Coefficient</b>, <b>Apex Score</b>, and <b>Raw elution profile</b>. "0" indicates that we don't use this correlation score and "1" indicates that we use this correlation score.
    * If you want to run Convolutional Neural Network (CNN) or Label Spreading (LS), you must set this parameter to "<b>`-s` `000000001`</b>". (* note that there are 9 characters in the string).
    * If you want to run EPPC with EPIC scores, then you can set this parameter to "<b>`-s`  `11101001`</b>". (* note that there are 8 characters in the string). In this example, it will use Mutual Information, Bayes Correlation, Euclidean Distance, Jaccard Score and Apex Score. To specify the correlation scores to use:

2. <b>`input_directory`</b>: This parameter stores the input directory where you store your elution profile file. It is recommended to use the abosulte path instead of relative path.


3.  (`-c` <b>`gold_standard_file_path`</b>) or (`--cluster` <b>`gold_standard_file_path`</b>): This parameter stores the path to the gold standard file that you curated.

4. <b>`output_directory`</b>: This parameter stores the path to the ouput directory. Make sure that you've already created the directory before running the command. It is recommended to use the abosulte path instead of relative path.

5. (`-o` <b>`output_filename_prefix`</b>) or (`--output_prefix` <b>`output_filename_prefix`</b>): You can specify a prefix name for all the output files. The default is "Out"

6. (`-M` <b>`training_method`</b>) or (`--classifier` <b>`training_method`</b>): This parameter specifies what kind of classifier that you use. Possible options include <b>`RF`</b>, <b>`CNN`</b>, <b>`LS`</b>. Note that <b>`RF`</b> must comes with selected EPIC scores like "<b>`-s`  `11101001`</b>" instead of raw elution profile ("<b>`-s` `000000001`</b>"). <b>`CNN`</b> and <b>`LS`</b> must come with raw elution profile ("<b>`-s` `000000001`</b>").

7. (`-n` <b>`number_of_cores`</b>) or (`--num_cores` <b>`number_of_cores`</b>): You need to specify the number of cores used to run EPPC, the default number is 1. Assume you want to use six cores to run EPIC, you can set "<b>`-n` `6`</b>"

8. `--LEARNING_SELECTION` <b>`learning method selection`</b>: This parameter specifies whether you want to use <b>supervised learning</b> or <b>semi-supervised learning</b>. If you want to run with <b>supervised learning</b>, then set "<b>`--LEARNING_SELECTION` `sl`</b>" (Your <b>`training_method`</b> can be <b>`RF`</b> or <b>`CNN`</b>); if you want to run with <b>semi-supervised learning</b>, then set "<b>`--LEARNING_SELECTION` `ssl`</b>" (Your <b>`training_method`</b> can be <b>`CNN`</b> or <b>`LS`</b>).

9. `--K_D_TRAIN` <b>`fold_or_direct_training`</b>: Set <b>`d`</b> to directly train the model; set <b>`k`</b> to run with k-fold training. (options: <b>`d`</b> and <b>`k`</b>; default: <b>`d`</b>)

10. `--FOLD_NUM` <b>`number_of_folds`</b>: If you set `--K_D_TRAIN` <b>`k`</b>, then this parameter stores how many folds you are going to evaluate your mode. Note that this parameter must be bigger than `2`. (default: <b>`5`</b>)

11. `--TRAIN_TEST_RATIO` <b>`testing_data_ratio`</b>: This parameter stores the ratio of testing data to all data. (default: <b>`0.3`</b>)

12. `--POS_NEG_RATIO` <b>`negative_PPIs_ratio`</b>: This parameter stores the ratio of negative PPIs to positive PPIs. (default: <b>`1`</b>)

13. `--NUM_EP` <b>`number_of_elution_profiles`</b>: This parameter stores the number of elution profiles inside each PPI. (default: <b>`2`</b>)

14. `--NUM_FRC` <b>`number_of_fractions`</b>:: This parameter stores the number of fractions in the elution profile file. (default: <b>`27`</b>)

## Supervised learning: EPIC scores + random Forest (EPIC scores + RF)

## Supervised Learning: raw inputs + convolutional neural networks (Raw + CNN)

## Semi-supervised learning: raw inputs + convolutional neural networks (Raw + CNN)

## Semi-supervised learning: raw inputs + label spreading (Raw + LS)

## Multiple 
