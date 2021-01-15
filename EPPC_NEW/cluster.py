import os
import csv

def clusters_prediction(ppis_file, cluster_file):
	clustering_CMD = "java -jar "+"./cluster_one-1.0.jar " + ppis_file + " > " + cluster_file
	os.system(clustering_CMD)
	# os.system(clustering_CMD)


class Clusters():
    def __init__(self, cl_path, file_type="tsv"):
        self.cl_path = cl_path
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
        cl_file = open(self.cl_path)
        cl_read = csv.reader(cl_file, delimiter=delimiter)
        for idx, line in enumerate(cl_read):
            self.add_complex(idx, line)
        cl_file.close()

    def set_prots_2_comps(self):
        for comp in self.complexes:
            for prot in self.complexes[comp]:
                if prot not in self.prots_2_comps:
                    self.prots_2_comps[prot] = set()
                self.prots_2_comps[prot].add(comp)
