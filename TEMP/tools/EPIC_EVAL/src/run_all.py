import os

# 1000, 100
n_estimator_ls = [1000]
# 20, 5, None
max_depth_ls = [20, 5, -1]
# 0.3
test_ratio_ls = [30]
# 1, 5
negative_ratio_ls = [1, 5]

eval_ratio_ls = [50]

for n_estimator in n_estimator_ls:
    for max_depth in max_depth_ls:
        for test_ratio in test_ratio_ls:
            for negative_ratio in negative_ratio_ls:
                for eval_ratio in eval_ratio_ls:
                    print("* n_estimator: ", n_estimator)
                    print("* max_depth: ", max_depth)
                    print("* test_ratio: ", test_ratio)
                    print("* negative_ratio: ", negative_ratio)
                    print("* eval_ratio: ", eval_ratio)

                    if (max_depth == -1):
                        max_depth_dir = "None"
                    else:
                        max_depth_dir = str(max_depth)

                    # exec_command = 'python ./main.py -s 11101001 /home/kuan-hao/EPPC/input/EPIC_test -c /home/kuan-hao/EPPC/input/Worm_reference_complexes.txt /home/kuan-hao/EPPC/output/EPIC_DATA/RF_EVAL/DIRECT/n_estimator_'+str(n_estimator)+'__max_depth_'+str(max_depth_dir)+'__test_ratio_'+str(test_ratio)+'__negative_ratio_'+str(negative_ratio)+'__eval_ratio_'+str(eval_ratio)+' -o Beads_A -M RF -n 10 -m EXP -f STRING --N_ESTIMATORS '+str(n_estimator)+' --MAX_DEPTH '+str(max_depth)+' --K_D_TRAIN d --FOLD_NUM 0 --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio)+' --EVAL_RATIO '+str(float(eval_ratio)/100)

                    exec_command = 'python ./src/main.py -s 11101001 /home/kuan-hao/EPPC/input/real_data/elution_profiles -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/SELF_DATA/RF_EVAL/DIRECT/n_estimator_'+str(n_estimator)+'__max_depth_'+str(max_depth_dir)+'__test_ratio_'+str(test_ratio)+'__negative_ratio_'+str(negative_ratio)+'__eval_ratio_'+str(eval_ratio)+' -o Intensity_H -M RF -n 10 -m EXP -f STRING --N_ESTIMATORS '+str(n_estimator)+' --MAX_DEPTH '+str(max_depth)+' --K_D_TRAIN d --FOLD_NUM 0 --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio)+' --EVAL_RATIO '+str(float(eval_ratio)/100)

                    print("exec_command: ", exec_command)
                    print("\n")
                    process = os.system(exec_command)
                # process = subprocess.Popen(['python', ])
                # print("process: ", process)
                # stdout, stderr = process.communicate()
                # print(stdout)
                # print(stderr)


# process = subprocess.Popen(['echo', 'More output'],
#                     stdout=subprocess.PIPE,
#                     stderr=subprocess.PIPE)
#
# stdout, stderr = process.communicate()
# print(stdout)
# print(stderr)