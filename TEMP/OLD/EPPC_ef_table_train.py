# Import classes from your brand new package
from EPPC_NEW import elution_profile_matrix as ep
from EPPC_NEW import gold_standard as gs
from EPPC_NEW import protein_protein_interaction as ppi
from EPPC_NEW import scores
from EPPC_NEW import classifier
import matplotlib.pyplot as plt

import itertools
import numpy as np
import pickle
import copy

ppi_pkl_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/EPIC_output/ppi.pkl"

with open(ppi_pkl_path, 'rb') as f:
    ppi = pickle.load(f)

print(len(ppi.positive))
print(len(ppi.negative))

## This is important for rebalancing the positive and negative ratio
rebalance_num = 5
ppi.rebalance(1/rebalance_num)

print(len(ppi.positive))
print(len(ppi.negative))

pcc_table = []
labels = []

for pos_ppi in ppi.positive:
    # print(np.array(pos_ppi.PCC_table_scores)[0].shape)
    diagonal_1 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(1), (0, 1), 'constant')
    diagonal_2 = np.array(pos_ppi.PCC_table_scores)[0].diagonal()
    diagonal_3 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(-1), (0, 1), 'constant')
    tmp_pcc_table = np.concatenate(([diagonal_1], [diagonal_2], [diagonal_3]), axis=0)
    pcc_table.append(tmp_pcc_table)
    labels.append(pos_ppi.label)

for pos_ppi in ppi.negative:
    diagonal_1 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(1), (0, 1), 'constant')
    diagonal_2 = np.array(pos_ppi.PCC_table_scores)[0].diagonal()
    diagonal_3 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(-1), (0, 1), 'constant')
    tmp_pcc_table = np.concatenate(([diagonal_1], [diagonal_2], [diagonal_3]), axis=0)
    pcc_table.append(tmp_pcc_table)
    labels.append(pos_ppi.label)

num_samples, row, frc_sz = np.array(pcc_table).shape
print("num_samples: ", num_samples)
print("row: ", row)
print("frc_sz: ", frc_sz)
labels = np.array(labels)

rnd_idx = np.arange(num_samples)
np.random.shuffle(rnd_idx)

# print(np.array(pcc_table).shape)
# print(rnd_idx)

test_sizes = [30]
# test_size = 30
for test_size in test_sizes:
    model = classifier.CNN_classifier(np.array(pcc_table)[rnd_idx], labels[rnd_idx], 200, "./output/09_16/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/pcc_table_diagonal/CNN/", "./output/09_16/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/pcc_table_diagonal/CNN/TEST/", "CNN_pcc_table_diagonal_model", test_size=test_size/100)

########################################
## Get whole positive & negative PPIS ##
########################################
    positive_PPI_name = []
    positive_eps = []
    positive_labels = []

    negative_PPI_name = []
    negative_eps = []
    negative_labels = []

    # unsure_PPI_name = []
    # unsure_eps = []
    # unsure_labels = []


    for pos_ppi in ppi.original_positive:
        positive_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))

        diagonal_1 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(1), (0, 1), 'constant')
        diagonal_2 = np.array(pos_ppi.PCC_table_scores)[0].diagonal()
        diagonal_3 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(-1), (0, 1), 'constant')
        tmp_pcc_table = np.concatenate(([diagonal_1], [diagonal_2], [diagonal_3]), axis=0)
        positive_eps.append(tmp_pcc_table)

        positive_labels.append(pos_ppi.label)

    for pos_ppi in ppi.original_negative:
        negative_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))

        diagonal_1 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(1), (0, 1), 'constant')
        diagonal_2 = np.array(pos_ppi.PCC_table_scores)[0].diagonal()
        diagonal_3 = np.pad(np.array(pos_ppi.PCC_table_scores)[0].diagonal(-1), (0, 1), 'constant')
        tmp_pcc_table = np.concatenate(([diagonal_1], [diagonal_2], [diagonal_3]), axis=0)
        negative_eps.append(tmp_pcc_table)

        negative_labels.append(pos_ppi.label)

    # for pos_ppi in ppi.original_unsure:
    #     unsure_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
    #     unsure_eps.append([pos_ppi.ef_a, pos_ppi.ef_b])
    #     unsure_labels.append(pos_ppi.label)


    positive_PPI_name = np.array(positive_PPI_name)
    positive_eps = np.array(positive_eps)
    positive_labels = np.array(positive_labels)

    negative_PPI_name = np.array(negative_PPI_name)
    negative_eps = np.array(negative_eps)
    negative_labels = np.array(negative_labels)

    # unsure_PPI_name = np.array(unsure_PPI_name)
    # unsure_eps = np.array(unsure_eps)
    # unsure_labels = np.array(unsure_labels)


    ###########################
    ## For CNN evaluation
    ###########################
    num_samples, num_scores, num_frc = np.array(positive_eps).shape
    positive_eps = positive_eps.reshape((num_samples, num_scores, num_frc, 1))
    positive_test_probas = model.predict(positive_eps)
    positive_test = copy.deepcopy(positive_test_probas)
    positive_test[positive_test >= 0.5] = 1
    positive_test[positive_test < 0.5] = 0

    # True positive : get prediction 1 [sensitivity]
    print("TP number: ", np.count_nonzero(positive_test == 1))
    print("TP ratio: ", np.count_nonzero(positive_test == 1) / len(positive_test))

    ################################
    ## Check these elution profiles
    ################################
    print("FN number: ", np.count_nonzero(positive_test == 0))
    print("FN ratio: ", np.count_nonzero(positive_test == 0) / len(positive_test))

    idx_1 = np.where(positive_test == 1)
    idx_0 = np.where(positive_test == 0)

    print(positive_eps.shape)
    print(positive_eps[idx_1[0]].shape)
    # os.makedirs("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TP/")
    # for i, eps in enumerate(positive_eps[idx_1[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     # if !os.exists("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TP/"):
    #     plt.savefig("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TP/"+str(i)+".png")
    #     plt.close()
    #     if i == 100:
    #         break

    print(positive_eps.shape)
    print(positive_eps[idx_0[0]].shape)
    # os.makedirs("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FN/")
    # for i, eps in enumerate(positive_eps[idx_0[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     plt.savefig("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FN/"+str(i)+".png")
    #     plt.close()
    #     if i == 100:
    #         break

    # positive_eps[positive_test == 1]

    # print(positive_eps == 1)
    # print(positive_eps[positive_test == 1])
    # print(positive_eps[positive_test == 0])




    num_samples, num_scores, num_frc = np.array(negative_eps).shape
    negative_eps = negative_eps.reshape((num_samples, num_scores, num_frc, 1))
    negative_test_probas = model.predict(negative_eps)
    negative_test = copy.deepcopy(negative_test_probas)
    negative_test[negative_test >= 0.5] = 1
    negative_test[negative_test < 0.5] = 0
    # False positive : get prediction 1
    # negative_test.count(1) / len(negative_test)
    print("FP number: ", np.count_nonzero(negative_test == 1))
    print("FP ratio: ", np.count_nonzero(negative_test == 1) / len(negative_test))

    ################################
    ## Check these elution profiles
    ################################
    print("TN number: ", np.count_nonzero(negative_test == 0))
    print("TN ratio: ", np.count_nonzero(negative_test == 0) / len(negative_test))



    idx_1 = np.where(negative_test == 1)
    idx_0 = np.where(negative_test == 0)

    print(negative_eps.shape)
    print(negative_eps[idx_1[0]].shape)
    # os.makedirs("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FP/")
    # for i, eps in enumerate(negative_eps[idx_1[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     plt.savefig("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FP/"+str(i)+".png")
    #     plt.close()
    #     if i == 100:
    #         break

    print(negative_eps.shape)
    print(negative_eps[idx_0[0]].shape)
    # os.makedirs("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TN/")
    # for i, eps in enumerate(negative_eps[idx_0[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     plt.savefig("./output/09_09/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TN/"+str(i)+".png")
    #     plt.close()
    #     if i == 100:
    #         break


    # if len(unsure_eps) != 0 :
    #     num_samples, num_scores, num_frc = np.array(unsure_eps).shape
    #     unsure_eps = unsure_eps.reshape((num_samples, num_scores, num_frc, 1))
    #     unsure_test_probas = model.predict(unsure_eps)
    #     unsure_test = copy.deepcopy(unsure_test_probas)
    #     unsure_test[unsure_test >= 0.5] = 1
    #     unsure_test[unsure_test < 0.5] = 0
    #     # Novel prediction : get prediction 1
    #     # unsure_test.count(1) / len(unsure_test)
    #     print("Novel number: ", np.count_nonzero(unsure_test == 1))
