import scipy.stats
from scipy.interpolate import interp1d
from scipy.spatial import distance
import math
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os

def intrpl_lin(prot_ef, intrpl_sz, frc_sz):
    x = np.linspace(0, frc_sz-1, num=frc_sz, endpoint=True)
    xnew = np.linspace(0, frc_sz-1, num=(intrpl_sz*(frc_sz-1)+frc_sz), endpoint=True)
    f = interp1d(x, prot_ef)
    return xnew, f(xnew)

def pearson_calculate(ef_a, ef_b):
    score = scipy.stats.pearsonr(ef_a, ef_b)[0]
    if math.isnan(score): return 0.0
    return score

def euclidean_dist_calculate(ef_a, ef_b):
    # print("ef_a: ", ef_a)
    # print("ef_b: ", ef_b)
    # print(all(ef_a == 0) and all(ef_b == 0))
    if all(ef_a == 0) and all(ef_b == 0):
        return 0
    else:
        return 1 - distance.euclidean(ef_a, ef_b)

class Pearson():
    def __init__(self, efm):
        self.name = "Pearson"
        self.efm = efm
        self.p2p_score = {}
        self.p2p_score_calculate()
        # self.score = self.pair_calculate(ef_a, ef_b)

    def p2p_score_calculate(self):
        # print(len(itertools.combinations(efm.ep_prot_list, 2)))
        for prot_pair in itertools.combinations(self.efm.ep_prot_list, 2):
            score = pearson_calculate(self.efm.get_elution_profile(prot_pair[0]), self.efm.get_elution_profile(prot_pair[1]))
            self.p2p_score[prot_pair] = score


class LocalPearson():
    def __init__(self, window_sz):
        self.name = "LocalPearson"
        self.window_sz = window_sz
        self.frac_scores = []
        # self.intrpl_epm = intrpl_epm
        # self.p2p_local_scores = {}
        # self.co_peaks = {}
        # self.pcc_score_cutoff = pcc_score_cutoff
        # self.peak_intst_cutoff = peak_intst_cutoff
        # self.p2p_score_calculate(1)
        # self.get_co_peak()

    def p2p_score_calculate(self, intrpl_ef_a, intrpl_ef_b, frc_sz, intrpl_sz):
        # intrpl_sz = self.intrpl_epm.intrpl_sz
        # frc_sz = self.intrpl_epm.frc_sz
        frac_scores = []
        window_sz = self.window_sz
        for idx in range(frc_sz):
            if (idx - window_sz) < 0:
                loc_intrpl_ef_a = intrpl_ef_a[0 : (idx+window_sz)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[0 : (idx+window_sz)*(intrpl_sz+1)+1]
            elif (idx + window_sz) > (frc_sz-1):
                loc_intrpl_ef_a = intrpl_ef_a[(idx-window_sz)*(intrpl_sz+1) : (frc_sz-1)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[(idx-window_sz)*(intrpl_sz+1) : (frc_sz-1)*(intrpl_sz+1)+1]
            else:
                loc_intrpl_ef_a = intrpl_ef_a[(idx-window_sz)*(intrpl_sz+1) : (idx+window_sz)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[(idx-window_sz)*(intrpl_sz+1) : (idx+window_sz)*(intrpl_sz+1)+1]
            loc_score = pearson_calculate(loc_intrpl_ef_a, loc_intrpl_ef_b)
            frac_scores.append(loc_score)
        self.frac_scores = frac_scores

class LocalPearsonTable():
    def __init__(self, window_sz):
        self.name = "LocalPearsonTable"
        self.window_sz = window_sz
        self.pcc_table = np.array([])
        # self.intrpl_epm = intrpl_epm
        # self.p2p_local_scores = {}
        # self.co_peaks = {}
        # self.pcc_score_cutoff = pcc_score_cutoff
        # self.peak_intst_cutoff = peak_intst_cutoff
        # self.p2p_score_calculate(1)
        # self.get_co_peak()

    def p2p_score_calculate(self, intrpl_ef_a, intrpl_ef_b, frc_sz, intrpl_sz):
        window_sz = self.window_sz
        loc_intrpl_ef_a_points = []
        loc_intrpl_ef_b_points = []

        exp_len = (intrpl_sz+1)*2+1
        for idx in range(frc_sz):
            if (idx - window_sz) < 0:
                loc_intrpl_ef_a = intrpl_ef_a[0 : (idx+window_sz)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[0 : (idx+window_sz)*(intrpl_sz+1)+1]
                if len(loc_intrpl_ef_a) < exp_len:
                    repeat_0s = np.repeat(0.0, exp_len-len(loc_intrpl_ef_a))
                    loc_intrpl_ef_a = np.concatenate([repeat_0s, np.array(loc_intrpl_ef_a)])
                    # print("len(loc_intrpl_ef_a): ", len(loc_intrpl_ef_a))
                    # print("loc_intrpl_ef_a: ", loc_intrpl_ef_a)
                if len(loc_intrpl_ef_b) < exp_len:
                    repeat_0s = np.repeat(0.0, exp_len-len(loc_intrpl_ef_b))
                    loc_intrpl_ef_b = np.concatenate([repeat_0s, np.array(loc_intrpl_ef_b)])
                    # print("len(loc_intrpl_ef_b): ", len(loc_intrpl_ef_b))
                    # print("loc_intrpl_ef_b: ", loc_intrpl_ef_b)
            elif (idx + window_sz) > (frc_sz-1):
                loc_intrpl_ef_a = intrpl_ef_a[(idx-window_sz)*(intrpl_sz+1) : (frc_sz-1)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[(idx-window_sz)*(intrpl_sz+1) : (frc_sz-1)*(intrpl_sz+1)+1]
                # print("len(loc_intrpl_ef_b): ", len(loc_intrpl_ef_b))
                # print("loc_intrpl_ef_b: ", loc_intrpl_ef_b)
                if len(loc_intrpl_ef_a) < exp_len:
                    repeat_0s = np.repeat(0.0, exp_len-len(loc_intrpl_ef_a))
                    loc_intrpl_ef_a = np.concatenate([np.array(loc_intrpl_ef_a), repeat_0s])
                    # print("len(loc_intrpl_ef_a): ", len(loc_intrpl_ef_a))
                    # print("loc_intrpl_ef_a: ", loc_intrpl_ef_a)
                if len(loc_intrpl_ef_b) < exp_len:
                    repeat_0s = np.repeat(0.0, exp_len-len(loc_intrpl_ef_b))
                    loc_intrpl_ef_b = np.concatenate([np.array(loc_intrpl_ef_b), repeat_0s])
                    # print("len(loc_intrpl_ef_b): ", len(loc_intrpl_ef_b))
                    # print("loc_intrpl_ef_b: ", loc_intrpl_ef_b)
            else:
                loc_intrpl_ef_a = intrpl_ef_a[(idx-window_sz)*(intrpl_sz+1) : (idx+window_sz)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[(idx-window_sz)*(intrpl_sz+1) : (idx+window_sz)*(intrpl_sz+1)+1]

            loc_intrpl_ef_a_points.append(loc_intrpl_ef_a)
            loc_intrpl_ef_b_points.append(loc_intrpl_ef_b)

        loc_intrpl_ef_a_points_len = len(loc_intrpl_ef_a_points)
        loc_intrpl_ef_b_points_len = len(loc_intrpl_ef_b_points)

        pcc_table = np.zeros((loc_intrpl_ef_a_points_len, loc_intrpl_ef_b_points_len))
        for i in range(loc_intrpl_ef_a_points_len):
            for j in range(loc_intrpl_ef_b_points_len):
                # print(i, j)
                loc_score = pearson_calculate(loc_intrpl_ef_a_points[i], loc_intrpl_ef_b_points[j])
                pcc_table[i][j] = loc_score
        # print("pcc_table: ", pcc_table)
        # print("pcc_table.shape: ", pcc_table.shape)
            # loc_score = pearson_calculate(loc_intrpl_ef_a, loc_intrpl_ef_b)
            # frac_scores.append(loc_score)
        self.pcc_table = pcc_table


class LocalEuclideanDist():
    def __init__(self, window_sz):
        self.name = "LocalEuclideanDist"
        self.window_sz = window_sz
        self.frac_scores = []

    def p2p_score_calculate(self, intrpl_ef_a, intrpl_ef_b, frc_sz, intrpl_sz):
        frac_scores = []
        window_sz = self.window_sz
        for idx in range(frc_sz):
            if (idx - window_sz) < 0:
                loc_intrpl_ef_a = intrpl_ef_a[0 : (idx+window_sz)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[0 : (idx+window_sz)*(intrpl_sz+1)+1]
            elif (idx + window_sz) > (frc_sz-1):
                loc_intrpl_ef_a = intrpl_ef_a[(idx-window_sz)*(intrpl_sz+1) : (frc_sz-1)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[(idx-window_sz)*(intrpl_sz+1) : (frc_sz-1)*(intrpl_sz+1)+1]
            else:
                loc_intrpl_ef_a = intrpl_ef_a[(idx-window_sz)*(intrpl_sz+1) : (idx+window_sz)*(intrpl_sz+1)+1]
                loc_intrpl_ef_b = intrpl_ef_b[(idx-window_sz)*(intrpl_sz+1) : (idx+window_sz)*(intrpl_sz+1)+1]
            loc_score = euclidean_dist_calculate(loc_intrpl_ef_a, loc_intrpl_ef_b)
            frac_scores.append(loc_score)
        self.frac_scores = frac_scores
