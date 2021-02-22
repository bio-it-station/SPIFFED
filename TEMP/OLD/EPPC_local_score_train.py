# Import classes from your brand new package
from EPPC_NEW import elution_profile_matrix as ep
from EPPC_NEW import gold_standard as gs
from EPPC_NEW import protein_protein_interaction as ppi
from EPPC_NEW import scores
from EPPC_NEW import classifier
from EPPC_NEW import cluster

import matplotlib.pyplot as plt

import itertools
import numpy as np
import pickle
import copy
import csv
import os

ppi_pkl_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/EPIC_output/ppi_add_unsure.pkl"

with open(ppi_pkl_path, 'rb') as f:
    ppi = pickle.load(f)

print(len(ppi.positive))
print(len(ppi.negative))

################################################
## Get the rebalanced positive & negative PPI ##
################################################
## This is important for rebalancing the positive and negative ratio
rebalance_num = 5
ppi.rebalance(1/rebalance_num)

print(len(ppi.positive))
print(len(ppi.negative))

scores = []
labels = []

for pos_ppi in ppi.positive:
    scores.append(pos_ppi.scores)
    labels.append(pos_ppi.label)

for pos_ppi in ppi.negative:
    scores.append(pos_ppi.scores)
    labels.append(pos_ppi.label)
    # labels = np.append(labels, pos_ppi.label)
    # scores = np.append(scores, pos_ppi.scores)

num_samples, num_scores, num_frc = np.array(scores).shape
print("num_samples: ", num_samples)
print("num_scores: ", num_scores)
print("num_frc: ", num_frc)
scores = np.array(scores)
labels = np.array(labels)
# print("scores: ", np.array(scores).shape)
# print("labels: ", np.array(labels).shape)

rnd_idx = np.arange(num_samples)
np.random.shuffle(rnd_idx)
# print("rnd_idx: ", len(rnd_idx))
# print("labels: ", len(labels))
# print("labels(rnd_idx): ", labels[rnd_idx])
#
# print("scores: ", len(scores))
# print("scores(rnd_idx): ", (scores[rnd_idx]).shape)

# test_sizes = [10, 15, 20, 25, 30]

#####################################
## Training model and benchmarking ##
#####################################
test_sizes = [30]
for test_size in test_sizes:
    # model = classifier.CNN_classifier(scores[rnd_idx], labels[rnd_idx], 200, "./output/09_16/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/local_scores/CNN/", "./output/09_16/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/local_scores/CNN/TEST", "CNN_scores_model", test_size = test_size/100)
    # svm = classifier.SVM_classifier(scores[rnd_idx], labels[rnd_idx], ppi.ratio, "./output/08_31/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/local_scores/SVM")
    clf = classifier.RF_classifier(scores[rnd_idx], labels[rnd_idx], ppi.ratio, "./output/09_16/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/local_scores/RF/", "./output/09_16/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/local_scores/RF/TEST", test_size = test_size/100)

    ###########################
    ## Original data
    ###########################

    # scores_tmp = scores[rnd_idx].reshape((num_samples, num_scores, num_frc, 1))
    # positive_test_probas = model.predict(scores_tmp)
    #
    # positive_test = copy.deepcopy(positive_test_probas)
    # positive_test[positive_test >= 0.5] = 1
    # positive_test[positive_test < 0.5] = 0
    #
    # # print((positive_test == 1))
    # # print((positive_test == scores_tmp))
    # # (positive_test == 1) && (positive_test == labels[rnd_idx])
    # positive_test = positive_test.reshape((num_samples,))
    # print("positive_test: ", positive_test)
    # print("labels[rnd_idx]: ", labels[rnd_idx])
    #
    # print(len(positive_test == 1 and labels[rnd_idx] == 1))
    # print(len(positive_test == 1 and labels[rnd_idx] == 0))
    # print(len(positive_test == 0 and labels[rnd_idx] == 1))
    # print(len(positive_test == 0 and labels[rnd_idx] == 0))
    #
    #
    # # True positive : get prediction 1 [sensitivity]
    # print("TP number: ", np.count_nonzero(positive_test == 1))
    # print("TP ratio: ", np.count_nonzero(positive_test == 1) / len(positive_test))
    #
    # print("FN number: ", np.count_nonzero(positive_test == 0))
    # print("FN ratio: ", np.count_nonzero(positive_test == 0) / len(positive_test))


    # ep_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/elution_profiles/beadsA.txt"
    # ep_data = ep.ElutionProfileMatrix(ep_path, file_type="tsv")
    # # print("ep_data: ", ep_data.ep_mat)
    # # intrpl_ep_data = ep.IntrplEPM(ep_path, 3, file_type="tsv")
    #
    # ########################################
    # ## Get whole positive & negative PPIS ##
    # ########################################
    # positive_PPI_name = []
    # positive_scores = []
    # positive_labels = []
    #
    # negative_PPI_name = []
    # negative_scores = []
    # negative_labels = []
    #
    # unsure_PPI_name = []
    # unsure_scores = []
    # unsure_labels = []
    #
    #
    # for pos_ppi in ppi.original_positive:
    #     positive_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
    #     positive_scores.append(pos_ppi.scores)
    #     positive_labels.append(pos_ppi.label)
    #
    # for pos_ppi in ppi.original_negative:
    #     negative_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
    #     negative_scores.append(pos_ppi.scores)
    #     negative_labels.append(pos_ppi.label)
    #
    # for pos_ppi in ppi.original_unsure:
    #     unsure_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
    #     unsure_scores.append(pos_ppi.scores)
    #     unsure_labels.append(pos_ppi.label)
    #
    #
    # positive_PPI_name = np.array(positive_PPI_name)
    # positive_scores = np.array(positive_scores)
    # positive_labels = np.array(positive_labels)
    #
    # negative_PPI_name = np.array(negative_PPI_name)
    # negative_scores = np.array(negative_scores)
    # negative_labels = np.array(negative_labels)
    #
    # unsure_PPI_name = np.array(unsure_PPI_name)
    # unsure_scores = np.array(unsure_scores)
    # unsure_labels = np.array(unsure_labels)
    #
    #
    # ###########################
    # ## For CNN evaluation
    # ###########################
    # num_samples, num_scores, num_frc = np.array(positive_scores).shape
    # positive_scores = positive_scores.reshape((num_samples, num_scores, num_frc, 1))
    # positive_test_probas = model.predict(positive_scores)
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
    # idx_1 = np.where(positive_test == 1)
    # idx_0 = np.where(positive_test == 0)
    #
    # # print(positive_scores.shape)
    # # print(positive_scores[idx_1[0]].shape)
    # # for i, eps in enumerate(positive_scores[idx_1[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/TP/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/TP/")
    # #     plt.plot(eps[0], 'x', color='blue');
    # #     # plt.plot(eps[1], 'o', color='black');
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/TP/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    # # for i, names in enumerate(positive_PPI_name[idx_1[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/TP/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/TP/")
    # #     plt.plot(ep_data.get_elution_profile(names[0]), '-', color='blue');
    # #     plt.plot(ep_data.get_elution_profile(names[1]), '-', color='black');
    # #     # if !os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/CNN/CNN/PCC_LOCAL_CHECK/TP/"):
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/TP/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    #
    #
    #
    # # print(positive_scores.shape)
    # # print(positive_scores[idx_0[0]].shape)
    # # for i, eps in enumerate(positive_scores[idx_0[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/FN/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/FN/")
    # #     plt.plot(eps[0], 'x', color='blue');
    # #     # plt.plot(eps[1], 'o', color='black');
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/FN/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    # # for i, names in enumerate(positive_PPI_name[idx_0[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/FN/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/FN/")
    # #     plt.plot(ep_data.get_elution_profile(names[0]), '-', color='blue');
    # #     plt.plot(ep_data.get_elution_profile(names[1]), '-', color='black');
    # #     # if !os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/CNN/CNN/PCC_LOCAL_CHECK/TP/"):
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/FN/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    #
    #
    # num_samples, num_scores, num_frc = np.array(negative_scores).shape
    # negative_scores = negative_scores.reshape((num_samples, num_scores, num_frc, 1))
    # negative_test_probas = model.predict(negative_scores)
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
    # idx_1 = np.where(negative_test == 1)
    # idx_0 = np.where(negative_test == 0)
    #
    # # print(negative_scores.shape)
    # # print(negative_scores[idx_1[0]].shape)
    # # for i, eps in enumerate(negative_scores[idx_1[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/FP/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/FP/")
    # #     plt.plot(eps[0], 'x', color='blue');
    # #     # plt.plot(eps[1], 'o', color='black');
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/FP/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    # # for i, names in enumerate(negative_PPI_name[idx_1[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/FP/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/FP/")
    # #     plt.plot(ep_data.get_elution_profile(names[0]), '-', color='blue');
    # #     plt.plot(ep_data.get_elution_profile(names[1]), '-', color='black');
    # #     # if !os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/CNN/CNN/PCC_LOCAL_CHECK/TP/"):
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/FP/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    # #
    # # print(negative_scores.shape)
    # # print(negative_scores[idx_0[0]].shape)
    # # for i, eps in enumerate(negative_scores[idx_0[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/TN/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/TN/")
    # #     plt.plot(eps[0], 'x', color='blue');
    # #     # plt.plot(eps[1], 'o', color='black');
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/PCC_LOCAL_CHECK/TN/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    # # for i, names in enumerate(negative_PPI_name[idx_0[0]]):
    # #     if not os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/TN/"):
    # #         os.makedirs("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/TN/")
    # #     plt.plot(ep_data.get_elution_profile(names[0]), '-', color='blue');
    # #     plt.plot(ep_data.get_elution_profile(names[1]), '-', color='black');
    # #     # if !os.path.exists("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/CNN/CNN/PCC_LOCAL_CHECK/TP/"):
    # #     plt.savefig("./output/09_16/1_" + str(rebalance_num) + "/" + str(test_sizes[0]) + "_test/local_scores/CNN/EF_CHECK/TN/"+str(i)+".png")
    # #     plt.close()
    # #     if i == 100:
    # #         break
    #
    #
    # if len(unsure_scores) != 0 :
    #     num_samples, num_scores, num_frc = np.array(unsure_scores).shape
    #     unsure_scores = unsure_scores.reshape((num_samples, num_scores, num_frc, 1))
    #     unsure_test_probas = model.predict(unsure_scores)
    #     unsure_test = copy.deepcopy(unsure_test_probas)
    #     unsure_test[unsure_test >= 0.5] = 1
    #     unsure_test[unsure_test < 0.5] = 0
    #     # Novel prediction : get prediction 1
    #     # unsure_test.count(1) / len(unsure_test)
    #     print("Novel number: ", np.count_nonzero(unsure_test == 1))
    # #     print("Novel ratio: ", np.count_nonzero(unsure_test == 1) / len(unsure_test))
    #
    # ###########################
    # ## For RF evaluation
    # ###########################
    # # num_samples, num_scores, num_frc = np.array(positive_scores).shape
    # # positive_scores = positive_scores.reshape((num_samples, num_scores*num_frc))
    # # positive_test = clf.predict(positive_scores)
    # # positive_test_probas = clf.predict_proba(positive_scores)[:,1]
    # # positive_test_probas_whole = clf.predict_proba(positive_scores)
    # # # True positive : get prediction 1 [sensitivity]
    # # print("TP number: ", np.count_nonzero(positive_test == 1))
    # # print("TP ratio: ", np.count_nonzero(positive_test == 1) / len(positive_test))
    # # print("FN number: ", np.count_nonzero(positive_test == 0))
    # # print("FN ratio: ", np.count_nonzero(positive_test == 0) / len(positive_test))
    # # print(len(positive_PPI_name[positive_test == 1]))
    # # positive_PPI_name_1 = positive_PPI_name[positive_test == 1]
    # # positive_test_probas_1 = positive_test_probas[positive_test == 1]
    # #
    # #
    # # num_samples, num_scores, num_frc = np.array(negative_scores).shape
    # # negative_scores = negative_scores.reshape((num_samples, num_scores*num_frc))
    # # negative_test = clf.predict(negative_scores)
    # # negative_test_probas = clf.predict_proba(negative_scores)[:,1]
    # # negative_test_probas_whole = clf.predict_proba(negative_scores)
    # # # True positive : get prediction 1 [sensitivity]
    # # print("FP number: ", np.count_nonzero(negative_test == 1))
    # # print("FP ratio: ", np.count_nonzero(negative_test == 1) / len(negative_test))
    # # print("TN number: ", np.count_nonzero(negative_test == 0))
    # # print("TN ratio: ", np.count_nonzero(negative_test == 0) / len(negative_test))
    # # print(len(negative_PPI_name[negative_test == 1]))
    # # negative_PPI_name_1 = negative_PPI_name[negative_test == 1]
    # # negative_test_probas_1 = negative_test_probas[negative_test == 1]
    # #
    # #
    # # if len(unsure_scores) != 0 :
    # #     num_samples, num_scores, num_frc = np.array(unsure_scores).shape
    # #     unsure_scores = unsure_scores.reshape((num_samples, num_scores*num_frc))
    # #     unsure_test = clf.predict(unsure_scores)
    # #     unsure_test_probas = clf.predict_proba(unsure_scores)[:,1]
    # #     unsure_test_probas_whole = clf.predict_proba(unsure_scores)
    # #     # Novel prediction : get prediction 1
    # #     # unsure_test.count(1) / len(unsure_test)
    # #     print("Novel number: ", np.count_nonzero(unsure_test == 1))
    # #     print("Novel ratio: ", np.count_nonzero(unsure_test == 1) / len(unsure_test))
    #
    #
    # # ###############################
    # # ## Writing PPI score into file
    # # ###############################
    # # output_ppi_file = "./output/ppi_interactions.txt"
    # # with open(output_ppi_file, 'w') as f:
    # #     writer = csv.writer(f, delimiter='\t')
    # #     # writer.writerow(['Protein_A', 'Protein_B', 'Score'])
    # #     for i in range(len(positive_PPI_name_1)):
    # #         writer.writerow(np.append(positive_PPI_name_1[i], str(positive_test_probas_1[i])))
    # #
    # # with open(output_ppi_file, 'w+') as f:
    # #     writer = csv.writer(f, delimiter='\t')
    # #     for i in range(len(negative_PPI_name_1)):
    # #         writer.writerow(np.append(negative_PPI_name_1[i], str(negative_test_probas_1[i])))
    # #
    # #
    # # ##################################
    # # ## Predicting clusters
    # # ##################################
    # # cluster.clusters_prediction("./output/ppi_interactions.txt", "./output/clusters_prediction.txt")
    #
    # ##################################
    # ## Evaluating predicted clusters
    # ##################################
    # # pred_clusters = GS.Clusters(False)
    # # pred_clusters.read_file("%s.clust.txt" % (output_dir))
    # # overlapped_complexes_with_reference = gs.get_complexes().get_overlapped_complexes_set(pred_clusters)
    # # print "# of complexes in reference dataset: " + str(len(overlapped_complexes_with_reference))
    # # #clust_scores, header = utils.clustering_evaluation(gs.complexes, pred_clusters, "", False)
    # # clust_scores, header, composite_score = utils.clustering_evaluation(gs.complexes, pred_clusters, "", False)
    # # outFH = open("%s.eval.txt" % (output_dir), "w")
    # # header = header.split("\t")
    # # clust_scores = clust_scores.split("\t")
    # # for i, head in enumerate(header):
    # # 	print "%s\t%s" % (head, clust_scores[i])
    # # 	print >> outFH, "%s\t%s" % (head, clust_scores[i])
    # # outFH.close()
