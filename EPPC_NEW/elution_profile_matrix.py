import csv
import numpy as np
from scipy.interpolate import interp1d
import os

class ElutionProfileMatrix():
    def __init__(self, ep_path, file_type='csv', frc_num_cutoff=2, max_frc_cutoff=1):
        self.ep_path = ep_path
        self.file_type = file_type
        self.frc_num_cutoff = frc_num_cutoff
        self.max_frc_cutoff = max_frc_cutoff
        self.ep_mat = []
        self.frc_sz = 0
        self.ep_prot_list = []
        self.ep_prot_dic = {}
        self.ep_removed_prot_list = []
        self.read_file()
         # = len(self.ep_mat[0])

    def read_file(self):
        paths = [os.path.join(self.ep_path,fn) for fn in next(os.walk(self.ep_path))[2]]
        ep_mats = {}
        ep_frc_szs = {}
        ep_prot_lists = {}
        ep_prot_dics = {}
        ep_removed_prot_lists = {}
        for epf in paths:
            if epf.rsplit(os.sep, 1)[-1].startswith("."): continue
            epf = epf.rstrip()
            print(epf)
            ep_file = open(epf)
            if self.file_type == "csv":
                delimiter = ','
            elif self.file_type == "tsv":
                delimiter = '\t'
            ep_read = csv.reader(ep_file, delimiter=delimiter)
            next(ep_read)
            idx_counter = 0
            lcl_ep_mat = {}
            lcl_ep_prot_list = []
            lcl_ep_prot_dic = {}
            lcl_ep_removed_prot_list = []
            for idx, line in enumerate(ep_read):
                prot_id = line[0]
                ep = np.array([float(i) for i in line[1:]])
                # print(len(ep))
                ep[np.isnan(ep)] = 0.0
                if (np.sum(ep > 0.0) > self.frc_num_cutoff) and (max(ep) > self.max_frc_cutoff):
                    lcl_ep_mat[prot_id] = ep
                    lcl_ep_prot_list.append(prot_id)
                    lcl_ep_prot_dic[prot_id] = idx_counter
                    idx_counter += 1
                else:
                    lcl_ep_removed_prot_list.append(prot_id)

            ep_mats[epf] = lcl_ep_mat
            ep_frc_szs[epf] = len(list(lcl_ep_mat.values())[0])
            ep_prot_lists[epf] = lcl_ep_prot_list
            ep_prot_dics[epf] = lcl_ep_prot_dic
            ep_removed_prot_lists[epf] = lcl_ep_removed_prot_list
            ep_file.close()
        print("ep_frc_szs: ", ep_frc_szs)
        # self.frc_sz = len(list(ep_mats[paths[0]].values())[0])
        # print("self.frc_sz: ", self.frc_sz)

        ########################
        ##### Now to merge #####
        ########################
        ep_prots = set()
        for k, v in ep_prot_dics.items():
            ep_prots = ep_prots | set(v.keys())
            # print(len(ep_prots))

        idx_counter = 0
        for prot in list(ep_prots):
            tmp_ep = []
            for epf in paths:
                if prot in ep_mats[epf]:
                    tmp_ep = np.concatenate(([tmp_ep], [ep_mats[epf][prot]]), axis=None)
                    # print(prot)
                else:
                    tmp_ep = np.concatenate(([tmp_ep], [np.repeat(0, ep_frc_szs[epf])]), axis=None)
                    # print("0")
            # print(tmp_ep)
            self.ep_mat.append(tmp_ep)
            self.ep_prot_list.append(prot)
            self.ep_prot_dic[prot] = idx_counter
            idx_counter += 1
        # print(self.ep_mat)
        # print(self.ep_prot_list)
        # print(self.ep_prot_dic)

        # for i in ep_mats:
        #     print(len(ep_prot_lists[i]))
        # for i in ep_prot_lists:
        #     print(len(ep_prot_lists[i]))
        # for i in ep_prot_dics:
        #     print(len(ep_prot_lists[i]))
        # for i in ep_removed_prot_lists:
        #     print(len(ep_prot_lists[i]))



        # print(ep_prot_lists)
# elutionDatas = []
# elutionProts = set([])
# for elutionFile in paths:
#     if elutionFile.rsplit(os.sep, 1)[-1].startswith("."): continue
#     elutionFile = elutionFile.rstrip()
#     elutionData = CS.ElutionData(elutionFile, frac_count=fc, max_frac_count=mfc)
#     if orthmap !="":
#         if orthmap != False:
#             mapper = GS.Inparanoid("", inparanoid_cutoff=1)
#             mapper.readTable(orthmap, direction=0)
#             elutionData.orthmap(mapper)
#     elutionDatas.append(elutionData)
#     elutionProts = elutionProts | set(elutionData.prot2Index.keys())
#     for score in scores:
#         score.init(elutionData)
# return elutionProts, elutionDatas



        print("Elution profile matrix '", self.ep_path,"' has been processed! ")
        print("       Total protein number: ", len(self.ep_removed_prot_list) + len(self.ep_prot_list))
        print("   Remaining protein number: ", len(self.ep_prot_list))
        # print("     Removed protein number: ", len(self.ep_removed_prot_list))
        # print(" Removed protein proportion: ", len(self.ep_removed_prot_list)/(len(self.ep_removed_prot_list) + len(self.ep_prot_list)))


    def has_prot(self, prot):
        if prot in self.ep_prot_list:
            return True
        else:
            return False

    def get_elution_profile(self, prot):
        if not self.has_prot(prot):
            return None
        else:
            return self.ep_mat[self.ep_prot_dic[prot]]


class IntrplEPM(ElutionProfileMatrix):
    def __init__(self, ep_path, intrpl_sz, intrpl_kind="linear", file_type='csv'):
        ElutionProfileMatrix.__init__(self, ep_path, file_type=file_type)
        self.intrpl_sz = intrpl_sz
        self.intrpl_kind = intrpl_kind
        self.intrpl_ep_mat = self.intrpl(intrpl_sz)

    def intrpl(self, intrpl_sz):
        intrpl_ep_mat = []
        x = np.linspace(0, self.frc_sz-1, num=self.frc_sz, endpoint=True)
        x_new = np.linspace(0, self.frc_sz-1, num=((intrpl_sz+1)*(self.frc_sz-1)+1), endpoint=True)
        for prot_ep in self.ep_mat:
            f = interp1d(x, prot_ep, kind=self.intrpl_kind)
            intrpl_ep_mat.append(f(x_new))
        return intrpl_ep_mat

    def get_intrpl_elution_profile(self, prot):
        if not self.has_prot(prot):
            return None
        else:
            return self.intrpl_ep_mat[self.ep_prot_dic[prot]]
