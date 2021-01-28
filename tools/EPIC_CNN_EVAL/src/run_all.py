import os

# 0.3
test_ratio_ls = [30]
# 1, 5
negative_ratio_ls = [1, 5]
K_D_TRAIN = 'd'
FOLD_NUM = 0
NUM_EP = 2
NUM_FRC = 60

eval_ratio_ls = [30, 50]

for test_ratio in test_ratio_ls:
    for negative_ratio in negative_ratio_ls:
        for eval_ratio in eval_ratio_ls:
            print("* test_ratio: ", test_ratio)
            print("* negative_ratio: ", negative_ratio)
            print("* eval_ratio: ", eval_ratio)
            print("* K_D_TRAIN: ", K_D_TRAIN)
            print("* FOLD_NUM: ", FOLD_NUM)
            print("* NUM_EP: ", NUM_EP)
            print("* NUM_FRC: ", NUM_FRC)


            exec_command = 'python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/EPIC_test -c /home/kuan-hao/EPPC/input/Worm_reference_complexes.txt /home/kuan-hao/EPPC/output/EPIC_DATA/CNN_EVAL/DIRECT/num_ep_'+str(NUM_EP)+'__num_frc_'+str(NUM_FRC)+'__test_ratio_'+str(test_ratio)+'__negative_ratio_'+str(negative_ratio)+'__eval_ratio_'+str(eval_ratio)+' -o Beads_A -M CNN -n 10 -m EXP -f STRING --K_D_TRAIN '+str(K_D_TRAIN)+' --FOLD_NUM '+str(FOLD_NUM)+' --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio) + ' --NUM_EP '+str(NUM_EP)+' --NUM_FRC '+str(NUM_FRC)+' --EVAL_RATIO '+str(float(eval_ratio)/100)
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