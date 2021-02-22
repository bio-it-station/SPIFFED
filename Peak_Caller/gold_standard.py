import os
import csv

# class GoldStandard():
#     def __init__(self, gs_path):
#         self.complexes = {}
#         # self.positive = ""
#         self.read_file(gs_path)

class Complexes():
    def __init__(self, gs_path, file_type="csv"):
        self.gs_path = gs_path
        self.file_type = file_type
        self.complexes = {}
        self.prots_2_comps = {}
        # self.positive = {}
        # self.negative = {}
        self.read_file()
        self.set_prots_2_comps()

    def add_complex(self, complex_name, prots):
        if complex_name in self.complexes.keys():
            self.complexes[complex_name] = set(self.complexes[complex_name] + prots)
        else:
            self.complexes[complex_name] = set(prots)

    def add_protein(self):
        if complex_name in self.complexes.keys():
            self.complexes[complex_name] = set(self.complexes[complex_name] + prots)
        else:
            self.complexes[complex_name] = set(prots)

    def read_file(self):
        if self.file_type == "csv":
            delimiter = ','
        elif self.file_type == "tsv":
            delimiter = '\t'
        gs_file = open(self.gs_path)
        gs_read = csv.reader(gs_file, delimiter=delimiter)
        for idx, line in enumerate(gs_read):
            self.add_complex(idx, line)
        gs_file.close()

    def set_prots_2_comps(self):
        for comp in self.complexes:
            for prot in self.complexes[comp]:
                if prot not in self.prots_2_comps:
                    self.prots_2_comps[prot] = set()
                self.prots_2_comps[prot].add(comp)
