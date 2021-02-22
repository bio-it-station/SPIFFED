# Import classes from your brand new package
from EPPC_NEW import elution_profile_matrix as ep
from EPPC_NEW import gold_standard as gs
from EPPC_NEW import protein_protein_interaction as ppi
from EPPC_NEW import scores
from EPPC_NEW import classifier
import matplotlib.pyplot as plt

import itertools
import numpy as np
import copy
import pickle
import os

ppi_pkl_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/EPIC_output/ppi_add_unsure.pkl"
ppi_pkl_save_new_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/EPIC_output/ppi_add_new.pkl"

with open(ppi_pkl_path, 'rb') as f:
    ppi = pickle.load(f)

print(len(ppi.positive))
print(len(ppi.negative))

## This is important for rebalancing the positive and negative ratio
rebalance_num = 3
test_size = 20

ppi.rebalance(1/rebalance_num)

print(len(ppi.positive))
print(len(ppi.negative))

with open(ppi_pkl_save_new_path, 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump(ppi, f)

efs = []
labels = []

for pos_ppi in ppi.positive:
    efs.append([pos_ppi.ef_a, pos_ppi.ef_b])
    labels.append(pos_ppi.label)

for pos_ppi in ppi.negative:
    efs.append([pos_ppi.ef_a, pos_ppi.ef_b])
    labels.append(pos_ppi.label)


# (4788, 2, 60)
num_samples, num_scores, num_frc = np.array(efs).shape
print("num_samples: ", num_samples)
print("num_scores: ", num_scores)
print("num_frc: ", num_frc)
efs = np.array(efs)
labels = np.array(labels)

rnd_idx = np.arange(num_samples)
np.random.shuffle(rnd_idx)


repeat_times = 20

for i in range(repeat_times):
    print(i)
    # for test_size in test_sizes:
    model = classifier.CNN_classifier(efs[rnd_idx], labels[rnd_idx], 200, "./output/09_13/"+str(i)+"/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/whole_ef/CNN/", "./output/09_13/"+str(i)+"/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/whole_ef/CNN/TEST/", "CNN_raw_ef_model", test_size = test_size/100)


    ########################################
    ## Get whole positive & negative PPIS ##
    ########################################
    positive_PPI_name = []
    positive_eps = []
    positive_labels = []

    negative_PPI_name = []
    negative_eps = []
    negative_labels = []

    unsure_PPI_name = []
    unsure_eps = []
    unsure_labels = []


    for pos_ppi in ppi.original_positive:
        positive_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
        positive_eps.append([pos_ppi.ef_a, pos_ppi.ef_b])
        positive_labels.append(pos_ppi.label)

    for pos_ppi in ppi.original_negative:
        negative_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
        negative_eps.append([pos_ppi.ef_a, pos_ppi.ef_b])
        negative_labels.append(pos_ppi.label)

    for pos_ppi in ppi.original_unsure:
        unsure_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
        unsure_eps.append([pos_ppi.ef_a, pos_ppi.ef_b])
        unsure_labels.append(pos_ppi.label)


    positive_PPI_name = np.array(positive_PPI_name)
    positive_eps = np.array(positive_eps)
    positive_labels = np.array(positive_labels)

    negative_PPI_name = np.array(negative_PPI_name)
    negative_eps = np.array(negative_eps)
    negative_labels = np.array(negative_labels)

    unsure_PPI_name = np.array(unsure_PPI_name)
    unsure_eps = np.array(unsure_eps)
    unsure_labels = np.array(unsure_labels)


    ###########################
    ## For CNN evaluation
    ###########################
    num_samples, num_scores, num_frc = np.array(positive_eps).shape
    positive_eps = positive_eps.reshape((num_samples, num_scores, num_frc, 1))
    positive_test_probas = model.predict(positive_eps)
    positive_test = copy.deepcopy(positive_test_probas)
    positive_test[positive_test >= 0.5] = 1
    positive_test[positive_test < 0.5] = 0
    # print("positive_test.shape[0]: ", positive_test.shape[0])
    tmp = positive_test.reshape((1, positive_test.shape[0]))
    print(tmp)
    positive_test_idx = np.where(tmp == 0)
    print("positive_test_idx: ", positive_test_idx)

    np.savetxt("./output/09_13/"+str(i)+"/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/whole_ef/CNN/FN_idx.txt", positive_test_idx[1].astype(int), newline=" ", fmt='% 4d')

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
    # os.makedirs("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TP/")
    # for i, eps in enumerate(positive_eps[idx_1[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     # if !os.exists("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TP/"):
    #     plt.savefig("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TP/"+str(i)+".png")
    #     plt.close()
    #     if i == 100:
    #         break
    #
    # print(positive_eps.shape)
    # print(positive_eps[idx_0[0]].shape)
    # os.makedirs("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FN/")
    # for i, eps in enumerate(positive_eps[idx_0[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     plt.savefig("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FN/"+str(i)+".png")
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
    # os.makedirs("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FP/")
    # for i, eps in enumerate(negative_eps[idx_1[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     plt.savefig("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/FP/"+str(i)+".png")
    #     plt.close()
    #     if i == 100:
    #         break
    #
    # print(negative_eps.shape)
    # print(negative_eps[idx_0[0]].shape)
    # os.makedirs("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TN/")
    # for i, eps in enumerate(negative_eps[idx_0[0]]):
    #     plt.plot(eps[0], '-', color='blue');
    #     plt.plot(eps[1], '-', color='black');
    #     plt.savefig("./output/09_13/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/whole_ef/CNN/EF_CHECK/TN/"+str(i)+".png")
    #     plt.close()
    #     if i == 100:
    #         break


    if len(unsure_eps) != 0 :
        num_samples, num_scores, num_frc = np.array(unsure_eps).shape
        unsure_eps = unsure_eps.reshape((num_samples, num_scores, num_frc, 1))
        unsure_test_probas = model.predict(unsure_eps)
        unsure_test = copy.deepcopy(unsure_test_probas)
        unsure_test[unsure_test >= 0.5] = 1
        unsure_test[unsure_test < 0.5] = 0
        # Novel prediction : get prediction 1
        # unsure_test.count(1) / len(unsure_test)
        print("Novel number: ", np.count_nonzero(unsure_test == 1))







    # ###########################
    # ## For CNN evaluation
    # ###########################
    # num_samples, num_scores, num_frc = np.array(positive_eps).shape
    # positive_eps = positive_eps.reshape((num_samples, num_scores, num_frc, 1))
    # positive_test_probas = model.predict(positive_eps)
    # positive_test = copy.deepcopy(positive_test_probas)
    # positive_test[positive_test >= 0.5] = 1
    # positive_test[positive_test < 0.5] = 0
    # # True positive : get prediction 1 [sensitivity]
    # print("TP number: ", np.count_nonzero(positive_test == 1))
    # print("TP ratio: ", np.count_nonzero(positive_test == 1) / len(positive_test))
    #
    # print("FN number: ", np.count_nonzero(positive_test == 0))
    # print("FN ratio: ", np.count_nonzero(positive_test == 0) / len(positive_test))
    #
    #
    # num_samples, num_scores, num_frc = np.array(negative_eps).shape
    # negative_eps = negative_eps.reshape((num_samples, num_scores, num_frc, 1))
    # negative_test_probas = model.predict(negative_eps)
    # negative_test = copy.deepcopy(negative_test_probas)
    # negative_test[negative_test >= 0.5] = 1
    # negative_test[negative_test < 0.5] = 0
    # # False positive : get prediction 1
    # # negative_test.count(1) / len(negative_test)
    # print("FP number: ", np.count_nonzero(negative_test == 1))
    # print("FP ratio: ", np.count_nonzero(negative_test == 1) / len(negative_test))
    #
    # print("TN number: ", np.count_nonzero(negative_test == 0))
    # print("TN ratio: ", np.count_nonzero(negative_test == 0) / len(negative_test))
    #
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
    #     print("Novel ratio: ", np.count_nonzero(unsure_test == 1) / len(unsure_test))

    # ###########################
    # ## For RF evaluation
    # ###########################
    # num_samples, num_scores, num_frc = np.array(positive_eps).shape
    # positive_eps = positive_eps.reshape((num_samples, num_scores*num_frc))
    # positive_test = clf.predict(positive_eps)
    # positive_test_probas = clf.predict_proba(positive_eps)[:,1]
    # positive_test_probas_whole = clf.predict_proba(positive_eps)
    # # True positive : get prediction 1 [sensitivity]
    # print("TP number: ", np.count_nonzero(positive_test == 1))
    # print("TP ratio: ", np.count_nonzero(positive_test == 1) / len(positive_test))
    # print(len(positive_PPI_name[positive_test == 1]))
    # positive_PPI_name_1 = positive_PPI_name[positive_test == 1]
    # positive_test_probas_1 = positive_test_probas[positive_test == 1]
    #
    #
    # num_samples, num_scores, num_frc = np.array(negative_eps).shape
    # negative_eps = negative_eps.reshape((num_samples, num_scores*num_frc))
    # negative_test = clf.predict(negative_eps)
    # negative_test_probas = clf.predict_proba(negative_eps)[:,1]
    # negative_test_probas_whole = clf.predict_proba(negative_eps)
    # # True positive : get prediction 1 [sensitivity]
    # print("FP number: ", np.count_nonzero(negative_test == 1))
    # print("FP ratio: ", np.count_nonzero(negative_test == 1) / len(negative_test))
    # print(len(negative_PPI_name[negative_test == 1]))
    # negative_PPI_name_1 = negative_PPI_name[negative_test == 1]
    # negative_test_probas_1 = negative_test_probas[negative_test == 1]
    #
    #
    # if len(unsure_eps) != 0 :
    #     num_samples, num_scores, num_frc = np.array(unsure_eps).shape
    #     unsure_eps = unsure_eps.reshape((num_samples, num_scores*num_frc))
    #     unsure_test = clf.predict(unsure_eps)
    #     unsure_test_probas = clf.predict_proba(unsure_eps)[:,1]
    #     unsure_test_probas_whole = clf.predict_proba(unsure_eps)
    #     # Novel prediction : get prediction 1
    #     # unsure_test.count(1) / len(unsure_test)
    #     print("Novel number: ", np.count_nonzero(unsure_test == 1))
    #     print("Novel ratio: ", np.count_nonzero(unsure_test == 1) / len(unsure_test))
    #
    #
    # ###############################
    # ## Writing PPI score into file
    # ###############################
    # output_ppi_file = "./output/ppi_interactions.txt"
    # with open(output_ppi_file, 'w') as f:
    #     writer = csv.writer(f, delimiter='\t')
    #     # writer.writerow(['Protein_A', 'Protein_B', 'Score'])
    #     for i in range(len(positive_PPI_name_1)):
    #         writer.writerow(np.append(positive_PPI_name_1[i], str(positive_test_probas_1[i])))
    #
    # with open(output_ppi_file, 'w+') as f:
    #     writer = csv.writer(f, delimiter='\t')
    #     for i in range(len(negative_PPI_name_1)):
    #         writer.writerow(np.append(negative_PPI_name_1[i], str(negative_test_probas_1[i])))
