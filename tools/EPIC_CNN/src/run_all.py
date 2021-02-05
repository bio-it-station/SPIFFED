import os


##########################################
## Parameters that you should change
##########################################
test_ratio_ls = [30, 70, 90]
negative_ratio_ls = [1, 5]
# 'k' or 'd'
K_D_TRAIN = 'd'
FOLD_NUM = 5
# 'RF_SL', 'CNN_SL', 'CNN_SSL', 'LS_SSL', or 'CNN_ENSEMBLE'
OUT_DIR_MTHD = 'CNN_SSL'
# 'sl': supervised learning; or 'ssl': semi-supervised learning
LEARNING_SELECTION = 'sl'
# 'SELF' or 'EPIC' or 'TEST'
data_source = 'SELF'
##########################################
##########################################

NUM_EP = 2

if OUT_DIR_MTHD == 'RF_SL':
    FEATURES = "11101001"
elif OUT_DIR_MTHD == 'CNN_SL' or OUT_DIR_MTHD == 'CNN_SSL' or OUT_DIR_MTHD == 'LS_SSL' or OUT_DIR_MTHD == 'CNN_ENSEMBLE':
    FEATURES = "000000001"

if data_source == 'SELF':
    EP_DIR = '/home/kuan-hao/EPPC/input/real_data/elution_profiles/'
    GD_FILE = '/home/kuan-hao/EPPC/input/real_data/gold_standard.tsv'
    OUT_DIR_SRC = 'SELF_DATA'
    NUM_FRC = 27
elif data_source == 'EPIC':
    EP_DIR = '/home/kuan-hao/EPPC/input/EPIC_test/'
    GD_FILE = '/home/kuan-hao/EPPC/input/Worm_reference_complexes.txt'
    OUT_DIR_SRC = 'EPIC_DATA'
    NUM_FRC = 60
elif data_source == 'TEST':
    EP_DIR = '/home/kuan-hao/EPPC/input/real_data/elution_profiles_test/'
    GD_FILE = '/home/kuan-hao/EPPC/input/real_data/gold_standard.tsv'
    OUT_DIR_SRC = 'TEST'
    NUM_FRC = 27

if K_D_TRAIN == 'd':
    OUT_DIR_TRN = 'DIRECT'
elif K_D_TRAIN == 'k':
    OUT_DIR_TRN = 'FOLDS'

for test_ratio in test_ratio_ls:
    for negative_ratio in negative_ratio_ls:
        print("* test_ratio: ", test_ratio)
        print("* negative_ratio: ", negative_ratio)
        print("* K_D_TRAIN: ", K_D_TRAIN)
        print("* FOLD_NUM: ", FOLD_NUM)
        print("* NUM_EP: ", NUM_EP)
        print("* NUM_FRC: ", NUM_FRC)

        # CNN with weights; EPIC data
        # exec_command = 'python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/EPIC_test -c /home/kuan-hao/EPPC/input/Worm_reference_complexes.txt /home/kuan-hao/EPPC/output/EPIC_DATA/CNN/FOLDS/num_ep_'+str(NUM_EP)+'__num_frc_'+str(NUM_FRC)+'__fold_number_'+str(FOLD_NUM)+'__negative_ratio_'+str(negative_ratio)+' -o Beads_A -M CNN -n 10 -m EXP -f STRING --K_D_TRAIN '+str(K_D_TRAIN)+' --FOLD_NUM '+str(FOLD_NUM)+' --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio) + ' --NUM_EP '+str(NUM_EP)+' --NUM_FRC '+str(NUM_FRC)

        # CNN with weights; SELF data
        # exec_command = 'python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/real_data/elution_profiles/ -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/SELF_DATA/CNN/FOLDS/num_ep_'+str(NUM_EP)+'__num_frc_'+str(NUM_FRC)+'__fold_number_'+str(FOLD_NUM)+'__negative_ratio_'+str(negative_ratio)+' -o Intensity_H -M CNN -n 10 -m EXP -f STRING --K_D_TRAIN '+str(K_D_TRAIN)+' --FOLD_NUM '+str(FOLD_NUM)+' --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio) + ' --NUM_EP '+str(NUM_EP)+' --NUM_FRC '+str(NUM_FRC)

        # CNN without weights; EPIC data
        exec_command = 'python ./main.py -s '+FEATURES+' '+EP_DIR+' -c '+GD_FILE+' /home/kuan-hao/EPPC/output/'+OUT_DIR_SRC+'/'+OUT_DIR_MTHD+'/'+OUT_DIR_TRN+'/num_ep_'+str(NUM_EP)+'__num_frc_'+str(NUM_FRC)+'__fold_number_'+str(FOLD_NUM)+'__negative_ratio_'+str(negative_ratio)+' -o Beads_A -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION '+str(LEARNING_SELECTION)+' --K_D_TRAIN '+str(K_D_TRAIN)+' --FOLD_NUM '+str(FOLD_NUM)+' --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio) + ' --NUM_EP '+str(NUM_EP)+' --NUM_FRC '+str(NUM_FRC)

        # CNN without weights; SELF data
        # exec_command = 'python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/real_data/elution_profiles/ -c /home/kuan-hao/EPPC/input/real_data/gold_standard.tsv /home/kuan-hao/EPPC/output/SELF_DATA/CNN_SSL/DIRECT_NO_WEIGHTS/2_num_ep_'+str(NUM_EP)+'__num_frc_'+str(NUM_FRC)+'__fold_number_'+str(FOLD_NUM)+'__negative_ratio_'+str(negative_ratio)+' -o Intensity_H -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION '+str(LEARNING_SELECTION)+' --K_D_TRAIN '+str(K_D_TRAIN)+' --FOLD_NUM '+str(FOLD_NUM)+' --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio) + ' --NUM_EP '+str(NUM_EP)+' --NUM_FRC '+str(NUM_FRC)

        print("exec_command: ", exec_command)
        print("\n")
        # process = oss.system(exec_command)
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
