import itertools
import timeit
import numpy as np
import random as rd
from .scores import LocalPearson, LocalEuclideanDist, LocalPearsonTable

class PPIS():
    def __init__(self, epm, gs, ratio=0.2):

        ## Original positive and negative PPI without filtering
        self.original_positive = set()
        self.original_negative = set()
        self.original_unsure = set()

        ## filtered positive and negative PPI
        self.positive = set()
        self.negative = set()
        self.unsure = set()

        # Positive / Negative
        self.ratio = ratio

        self.init(epm, gs)
        self.rebalance(self.ratio)


    def init(self, epm, gs):
        start = timeit.default_timer()
        print("gs.prots_2_comps.keys(): ", len(gs.prots_2_comps.keys()))
        counter = 0
        for prot_pair in itertools.combinations(gs.prots_2_comps.keys(), 2):
            counter += 1
            ef_a = epm.get_elution_profile(prot_pair[0])
            ef_b = epm.get_elution_profile(prot_pair[1])

            # ef_a_intrpl = epm.get_intrpl_elution_profile(prot_pair[0])
            # ef_b_intrpl = epm.get_intrpl_elution_profile(prot_pair[1])
            ef_a_intrpl = []
            ef_b_intrpl = []
            # If the elution profile of the protein cannot be found
            if ef_a is None or ef_b is None: continue

            # If the two elution profile fractions has no intersection
            if np.all(~((ef_a!=0.0) & (ef_b!=0.0))):
                # print("     No interactions between fractions")
                # print("         ", ef_a)
                # print("         ", ef_b)
                continue

            ppi = PPI(prot_pair[0], prot_pair[1], ef_a, ef_b, ef_a_intrpl, ef_b_intrpl, epm, gs)
            if ppi.label == 1:
                self.original_positive.add(ppi)
            elif ppi.label == 0:
                self.original_negative.add(ppi)
            elif ppi.label is None:
                self.original_unsure.add(ppi)
                self.unsure.add(ppi)

        print("Counter: ", counter)
        stop = timeit.default_timer()
        print('Time: ', stop - start)

    # def get_ppi_lable(self, ppi_id):
    #     return labels_dic[ppi_id]
    #
    # def get_ppi_score(self, ppi_id):
    #     pass scores_dic[ppi_id]

    # def set_pos_neg(self, gs):
    #     for prot_pair in itertools.combinations(gs.prots_2_comps.keys, 2):
    #         if bool(gs.prots_2_comps[prot_pair[0]] & gs.prots_2_comps[prot_pair[1]]):
    #             self.positive.add(prot_pair)
    #         else:
    #             self.negative.add(prot_pair)

    def rebalance(self, ratio):
        if len(self.original_positive)/self.ratio > len(self.original_negative):
            self.positive = set(rd.sample(self.original_positive, int(len(self.original_negative)*ratio)))
            self.negative = self.original_negative
        else:
            self.positive = self.original_positive
            self.negative = set(rd.sample(self.original_negative, int(len(self.original_positive)/ratio)))

class PPI():
    def __init__(self, prot_a, prot_b, ef_a, ef_b, ef_a_intrpl, ef_b_intrpl, epm, gs):
        self.prot_a = prot_a
        self.prot_b = prot_b
        self.ef_a = ef_a
        self.ef_b = ef_b
        self.ef_a_intrpl = ef_a_intrpl
        self.ef_b_intrpl = ef_b_intrpl
        self.label = None
        self.scores = []
        self.scores_name = []
        self.PCC_table_scores = []
        self.init(prot_a, prot_b, ef_a, ef_b, ef_a_intrpl, ef_b_intrpl, epm, gs)

    def init(self, prot_a, prot_b, ef_a, ef_b, ef_a_intrpl, ef_b_intrpl, epm, gs):
        ## Pearson Correlation Table
        # lcl_pcc = LocalPearsonTable(window_sz=1)
        # lcl_pcc.p2p_score_calculate(ef_a_intrpl, ef_b_intrpl, epm.frc_sz, epm.intrpl_sz)
        # self.PCC_table_scores.append(lcl_pcc.pcc_table)
        #
        # ## Local Pearson Correlation
        # lcl_pcc = LocalPearson(window_sz=1)
        # lcl_pcc.p2p_score_calculate(ef_a_intrpl, ef_b_intrpl, epm.frc_sz, epm.intrpl_sz)
        # self.scores_name.append('LocalPearson')
        # self.scores.append(lcl_pcc.frac_scores)
        #
        # ## Local Euclidean Distance
        # lcl_euc = LocalEuclideanDist(window_sz=1)
        # lcl_euc.p2p_score_calculate(ef_a_intrpl, ef_b_intrpl, epm.frc_sz, epm.intrpl_sz)
        # self.scores_name.append('LocalEuclideanDist')
        # self.scores.append(lcl_euc.frac_scores)

        if (prot_a in gs.prots_2_comps) and (prot_b in gs.prots_2_comps):
            # Both proteins are in the gold standard
            if bool(gs.prots_2_comps[prot_a] & gs.prots_2_comps[prot_b]):
                # Proteins are in the same complex
                self.label = 1
            else:
                # Proteins are not in the same complex
                self.label = 0
        else:
            # Maybe one or both proteins are not in the gold standard
            self.label = None
