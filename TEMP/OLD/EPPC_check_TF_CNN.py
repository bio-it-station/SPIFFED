# Import classes from your brand new package
from EPPC_NEW import elution_profile_matrix as ep
from EPPC_NEW import gold_standard as gs
from EPPC_NEW import protein_protein_interaction as ppi
from EPPC_NEW import scores
from EPPC_NEW import classifier
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
import csv
import pickle
from dtw import dtw

ppi_pkl_save_new_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/EPIC_output/ppi_add_new.pkl"

with open(ppi_pkl_save_new_path, 'rb') as f:
    ppi = pickle.load(f)


print("Length: ", len(ppi.positive))
#
#
# ########################################
# ## Get whole positive & negative PPIS ##
# ########################################
# positive_PPI_name = []
# positive_eps = []
# positive_labels = []
#
# negative_PPI_name = []
# negative_eps = []
# negative_labels = []
#
# unsure_PPI_name = []
# unsure_eps = []
# unsure_labels = []
#
# for pos_ppi in ppi.positive:
#     positive_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
#     positive_eps.append([pos_ppi.ef_a, pos_ppi.ef_b])
#     positive_labels.append(pos_ppi.label)
#
# for pos_ppi in ppi.negative:
#     negative_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
#     negative_eps.append([pos_ppi.ef_a, pos_ppi.ef_b])
#     negative_labels.append(pos_ppi.label)
#
# for pos_ppi in ppi.original_unsure:
#     unsure_PPI_name.append((pos_ppi.prot_a, pos_ppi.prot_b))
#     unsure_eps.append([pos_ppi.ef_a, pos_ppi.ef_b])
#     unsure_labels.append(pos_ppi.label)
#
#
# positive_PPI_name = np.array(positive_PPI_name)
# positive_eps = np.array(positive_eps)
# positive_labels = np.array(positive_labels)
#
# negative_PPI_name = np.array(negative_PPI_name)
# negative_eps = np.array(negative_eps)
# negative_labels = np.array(negative_labels)
#
# unsure_PPI_name = np.array(unsure_PPI_name)
# unsure_eps = np.array(unsure_eps)
# unsure_labels = np.array(unsure_labels)
#
#
# ###########################
# ## For CNN evaluation
# ###########################
# num_samples, num_scores, num_frc = np.array(positive_eps).shape
# positive_eps = positive_eps.reshape((num_samples, num_scores, num_frc, 1))
# positive_test_probas = model.predict(positive_eps)
# positive_test = copy.deepcopy(positive_test_probas)
#
# # True positive : get prediction 1 [sensitivity]
# print("TP number: ", np.count_nonzero(positive_test == 1))
# print("TP ratio: ", np.count_nonzero(positive_test == 1) / len(positive_test))
#
# ################################
# ## Check these elution profiles
# ################################
# print("FN number: ", np.count_nonzero(positive_test == 0))
# print("FN ratio: ", np.count_nonzero(positive_test == 0) / len(positive_test))
#
# idx_1 = np.where(positive_test == 1)
# idx_0 = np.where(positive_test == 0)
#
# print(positive_eps.shape)
# print(positive_eps[idx_1[0]].shape)











rebalance_num = 3
test_size = 20
repeat_times = 20

global_counter = np.array([])
global_length = 0

# ppi.positive

for i in range(repeat_times):
    # for test_size in test_sizes:
    rd = np.loadtxt("./output/09_13/"+str(i)+"/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/whole_ef/CNN/FN_idx.txt")
    global_length = global_length + len(rd)
    # print(len(rd))
    global_counter = np.concatenate((global_counter, rd), axis=None)
    # with open("./output/09_13/"+str(i)+"/CNN_TUNE/1_" + str(rebalance_num) + "/" + str(test_size) + "_test/whole_ef/CNN/FN_idx.txt") as fd:
    #     rd = csv.reader(fd, delimiter="\t", quotechar='"')
    #     for row in rd:
    #         print(row)

counter_list = []
for i in range(623):
    count_num = np.count_nonzero(global_counter == i)
    # print(i, ": ", count_num)
    counter_list.append(count_num)

# print("counter_list: ", counter_list)
# print("counter_list length: ", len(counter_list))

counter_list_idx = np.where(np.array(counter_list) >= 18)
print(counter_list_idx[0])

if not os.path.exists("./output/09_13/ef_FP/"):
    os.mkdir("./output/09_13/ef_FP/")

ppi_positive_list = np.array(list(ppi.positive))
print(ppi_positive_list[[1, 2]])

# print("ppi_positive_list[counter_list_idx[0]]: ", ppi_positive_list[counter_list_idx[0]])
counter = 1
for i in ppi_positive_list[counter_list_idx[0]]:
    ef_a = i.ef_a
    ef_b = i.ef_b
    plt.plot(ef_a, '-', color='blue');
    plt.plot(ef_b, '-', color='black');
    plt.savefig("./output/09_13/ef_FP/" + str(counter) + ".png")
    plt.close()

    manhattan_distance = lambda ef_a, ef_b: np.abs(ef_a - ef_b)
    d, cost_matrix, acc_cost_matrix, path = dtw(ef_a, ef_b, dist=manhattan_distance)
    plt.imshow(acc_cost_matrix.T, origin='lower', cmap='gray', interpolation='nearest')
    plt.plot(path[0], path[1], 'w')
    plt.savefig("./output/09_13/ef_FP/" + str(counter) + "_dtw.png")
    plt.close()

    counter += 1
    # print("ef_a: ", ef_a)
    # print("ef_b: ", ef_b)


# print(ppi_positive_list[counter_list_idx[0]])

plt.bar(range(623), counter_list)
plt.savefig("./output/09_13/bar_plot.png", dpi=400)
plt.close()
