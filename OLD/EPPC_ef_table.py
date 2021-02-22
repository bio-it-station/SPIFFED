# Import classes from your brand new package
from EPPC_NEW import elution_profile_matrix as ep
from EPPC_NEW import gold_standard as gs
from EPPC_NEW import protein_protein_interaction as ppi

# from EPPC_NEW import protein_protein_interaction as ppi
# from EPPC_NEW import scores
# import matplotlib.pyplot as plt

import itertools
import numpy as np
import pickle

ep_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/elution_profiles/beadsA.txt"
gs_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/Worm_reference_complexes.txt"
ep_data = ep.ElutionProfileMatrix(ep_path, file_type="tsv")
# print("ep_data: ", ep_data.ep_mat)

intrpl_ep_data = ep.IntrplEPM(ep_path, 3, file_type="tsv")
gs_data= gs.Complexes(gs_path, file_type="tsv")
print(len(gs_data.complexes))
# print(gs_data.prots_2_comps)
print(len(gs_data.prots_2_comps))

ppi = ppi.PPIS(intrpl_ep_data, gs_data, ratio=0.2)
print(len(ppi.positive))
print(len(ppi.negative))

ppi_pkl_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/EPIC_output/ppi_table.pkl"

with open(ppi_pkl_path, 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump(ppi, f)
# print("ppi.labels: ", np.array(ppi.labels).shape)
# print("ppi.scores: ", np.array(ppi.scores).shape)
#
# positive = 0
# negative = 0
# for key, value in ppi.labels_dic.items():
#     if value:
#         positive += 1
#     elif not value:
#         negative += 1
# print("Positive: ", positive)
# print("Negative: ", negative)
#
# ppi_pkl_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/EPPC_NEW_dev_env/EPIC_input_test/EPIC_output/ppi.pkl"
#
# with open(ppi_pkl_path, 'wb') as f:  # Python 3: open(..., 'wb')
#     pickle.dump(ppi, f)
