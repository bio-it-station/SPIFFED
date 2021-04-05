import os

# python ./main.py -s 000000001 /home/kuan-hao/EPPC/input/OUR_DATA/intensity_high/ -c /home/kuan-hao/EPPC/input/OUR_DATA/gold_standard.tsv /home/kuan-hao/EPPC/output/TEST/ -o TEST -M CNN -n 6 -m EXP -f STRING --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 1 --NUM_EP 2 --NUM_FRC 27

##########################################
## Parameters that you should change
##########################################
## This is for random forest

test_ratio_ls = [30]
negative_ratio_ls = [1, 2, 5]
# 'k' or 'd'
K_D_TRAIN = 'k'
FOLD_NUM = 5
# 'RF_SL', 'CNN_SL', 'CNN_SSL', 'LS_SSL', or 'CNN_ENSEMBLE'
OUT_DIR_MTHD = 'CNN_SL_EVAL'
# 'sl': supervised learning; or 'ssl': semi-supervised learning
LEARNING_SELECTION = 'sl'
# 'SELF' or 'EPIC' or 'TEST'
data_source = 'EPIC'
##########################################
##########################################

NUM_EP = 2

if K_D_TRAIN == 'd':
    FOLD_NUM = '0'

if OUT_DIR_MTHD == 'RF_SL':
    FEATURES = "11101001"
    CNN_ENSEMBLE = 0

elif OUT_DIR_MTHD == 'CNN_SL_EVAL':
    FEATURES = "000000001"
    if OUT_DIR_MTHD == 'CNN_ENSEMBLE':
        CNN_ENSEMBLE = 1
    else:
        CNN_ENSEMBLE = 0

if data_source == 'SELF':
    EP_DIR = '/home/kuan-hao/EPPC/input/OUR_DATA/'
    GD_FILE = '/home/kuan-hao/EPPC/input/OUR_DATA/gold_standard.tsv'
    OUT_DIR_SRC = 'SELF_DATA'
    NUM_FRC = 27
    OUT_PREFIXS = ['intensity_high', 'intensity_medium', 'intensity_low']

elif data_source == 'EPIC':
    EP_DIR = '/home/kuan-hao/EPPC/input/EPIC_DATA/'
    GD_FILE = '/home/kuan-hao/EPPC/input/EPIC_DATA/Worm_reference_complexes.txt'
    OUT_DIR_SRC = 'EPIC_DATA'
    NUM_FRC = 60
    OUT_PREFIXS = ['beadsA', 'beadsB', 'beadsBNF', 'beadsC', 'beadsCRF', 'beadsL']

elif data_source == 'TEST':
    EP_DIR = '/home/kuan-hao/EPPC/input/OUR_DATA/elution_profiles_test/'
    GD_FILE = '/home/kuan-hao/EPPC/input/OUR_DATA/gold_standard.tsv'
    OUT_DIR_SRC = 'TEST'
    NUM_FRC = 27
    OUT_PREFIX = 'TEST'

if K_D_TRAIN == 'd':
    OUT_DIR_TRN = 'DIRECT'
elif K_D_TRAIN == 'k':
    OUT_DIR_TRN = 'FOLDS'

if OUT_DIR_MTHD == 'CNN_SL_EVAL':
    for test_ratio in test_ratio_ls:
        for negative_ratio in negative_ratio_ls:
            for OUT_PREFIX in OUT_PREFIXS:
                print("* test_ratio: ", test_ratio)
                print("* negative_ratio: ", negative_ratio)
                print("* K_D_TRAIN: ", K_D_TRAIN)
                print("* FOLD_NUM: ", FOLD_NUM)
                print("* NUM_EP: ", NUM_EP)
                print("* NUM_FRC: ", NUM_FRC)

                # CNN without weights; EPIC data
                exec_command = 'python ../../main.py -s '+FEATURES+' '+EP_DIR+OUT_PREFIX+' -c '+GD_FILE+' /home/kuan-hao/EPPC/output/'+OUT_DIR_SRC+'/'+OUT_PREFIX+'/'+OUT_DIR_MTHD+'/'+OUT_DIR_TRN+'/'+OUT_PREFIX+'__K_D__'+K_D_TRAIN+'__'+OUT_DIR_MTHD+'__num_ep_'+str(NUM_EP)+'__num_frc_'+str(NUM_FRC)+'__fold_number_'+str(FOLD_NUM)+'__negative_ratio_'+str(negative_ratio)+'__test_ratio_'+str(test_ratio)+' -o '+OUT_PREFIX+' -M CNN -n 10 -m EXP -f STRING --K_D_TRAIN '+str(K_D_TRAIN)+' --FOLD_NUM '+str(FOLD_NUM)+' --TRAIN_TEST_RATIO '+str(float(test_ratio)/100)+' --POS_NEG_RATIO '+str(negative_ratio) + ' --NUM_EP '+str(NUM_EP)+' --NUM_FRC '+str(NUM_FRC)

                print("exec_command: ", exec_command)
                print("\n")
                process = os.system(exec_command)
                # process = subprocess.Popen(['python', ])
                # print("process: ", process)
                # stdout, stderr = process.communicate()
                # print(stdout)
                # print(stderr)
