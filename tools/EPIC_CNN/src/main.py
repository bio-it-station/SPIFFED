from __future__ import division

import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
import sys
import copy
import os
import argparse
import pickle
import numpy as np

import warnings
warnings.filterwarnings('ignore')

# def Goldstandard_from_cluster_File(gsF, foundprots=""):
# 	clusters = GS.Clusters(need_to_be_mapped=False)
# 	clusters.read_file(gsF)
# 	if foundprots != "": clusters.remove_proteins(foundprots)
# 	gs = GS.Goldstandard_from_Complexes("All")
# 	gs.complexes = clusters
# 	gs.make_pos_neg_ppis()
# 	return gs
#
#
# def Goldstandard_from_PPI_File(gsF, foundprots=""):
# 	out = GS.Goldstandard_from_Complexes("gs")
# 	gsFH = open(gsF)
# 	for line in gsFH:
# 		line = line.rstrip()
# 		ida, idb, class_label = line.split("\t")[0:3]
# 		if foundprots !="" and (ida not in foundprots or idb not in foundprots): continue
# 		edge = "\t".join(sorted([ida, idb]))
# 		if class_label == "positive":
# 			out.positive.add(edge)
# 		else:
# 			out.negative.add(edge)
# 	gsFH.close()
# 	return out

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--feature_selection", type = str, help="Select which features to use. This is an 8 position long array of 0 and 1, where each position determines which co-elution feature to use. Features sorted by position are: MI, Bayes, Euclidean, WCC, Jaccard, PCCN, PCC, and Apex.  Each default=11101001", default="11101001")
	parser.add_argument("input_dir",  type = str, help="Directory containing the elution files for each experiment")

	parser.add_argument("-t", "--taxid", type = str, help="TAXID to automatically download reference from GO,CORUM,INtACT",
						default="")
	parser.add_argument("-c", "--cluster", type = str, help="Path to file containing protein clsuter reference",
						default="")
	parser.add_argument("-p", "--ppi", type = str, help="path to ppi File",
						default="")

	parser.add_argument("output_dir", type = str,help="Directory containing the output files")
	parser.add_argument("-o", "--output_prefix", type = str,help="Prefix name for all output Files", default="Out")

	parser.add_argument("-M", "--classifier", type = str,help="Select which classifier to use. Values: RF SVM, default RF",
						default="RF")
	parser.add_argument("-n", "--num_cores", type = int,help="Number of cores to be used, default 1",
						default=1)

	parser.add_argument("-m", "--mode", type=str,
						help="Run EPIC with experimental, functional, or both evidences. Values: EXP, FA, COMB, default: EXP  ",
						default="EXP")
	parser.add_argument("-f", "--fun_anno_source", type = str,help="Where to get functional annotaiton from. Values: STRING or GM or FILE, default= GM",
						default="GM")
	parser.add_argument("-F", "--fun_anno_file", type=str,
						help="Path to File containing functional annotation. This flag needs to be set when using FILE as fun_anno_source.",
						)
	parser.add_argument("-r", "--co_elution_cutoff", type = float,help="Co-elution score cutoff. default 0.5",
						default=0.5)
	parser.add_argument("-R", "--classifier_cutoff", type = float,help="Classifier confidence valye cutoff. default = 0.5",
						default=0.5)
	parser.add_argument("-e", "--elution_max_count", type = int,help="Removies protein that have a maximal peptide count less than the given value. default = 1",
						default=1)
	parser.add_argument("-E", "--frac_count", type = int,help="Number of fracrions a protein needs to be measured in. default = 2",
						default=2)

	parser.add_argument("-P", "--precalcualted_score_file", type = str,help="Path to precalulated scorefile to read scores from for faster rerunning of EPIC. default = None",
						default="NONE")

	########################
	## New added argument ##
	########################
	# k-fold or direct training:
	parser.add_argument("--K_D_TRAIN", type = str, help="Values: d, k", default="d")
	# fold_num
	parser.add_argument("--FOLD_NUM", type = int, help="default = 5", default=5)
	# train_test_ratio
	parser.add_argument("--TRAIN_TEST_RATIO", type = float, help="default = 0.3", default=0.3)
	# pos_neg_ratio
	parser.add_argument("--POS_NEG_RATIO", type = int, help="default = 1", default=1)
	# num_ep
	parser.add_argument("--NUM_EP", type = int, help="default = 2", default=2)
	# num_frc
	parser.add_argument("--NUM_FRC", type = int, help="default = 27", default=27)
	############################
	## New added argument END ##
	############################

	args = parser.parse_args()

	args.mode = args.mode.upper()
	args.fun_anno_source = args.fun_anno_source.upper()

	##########################################
	## This part is "PASE INPUT"
	##########################################
	# K_D_TRAIN = raw_input("Enter k-fold or direct training (k/d):")
	# if K_D_TRAIN == "k":
	# 	print "You are using k-fold training!"
	# 	FOLD_NUM = input("Enter fold number:")
	# 	print("Fold number is: " + str(FOLD_NUM))
	# 	TRAIN_TEST_RATIO = None
	# elif K_D_TRAIN == "d":
	# 	print "You are using direct training!"
	# 	TRAIN_TEST_RATIO = input("Enter test/all_data ratio:")
	# 	print("Your TRAIN_TEST_RATIO ratio is: " + str(TRAIN_TEST_RATIO))
	# 	FOLD_NUM = None
	# else:
	# 	print "Please input 'k' or 'd'."
	# 	sys.exit()
	# NUM_EP = 2
	# NUM_FRC = input("Enter number of fractions:")
	# POS_NEG_RATIO = input("Enter negative/positive ratio:")
	# num_ep = 2
	# num_frc = 27
	K_D_TRAIN = args.K_D_TRAIN
	FOLD_NUM = args.FOLD_NUM
	if FOLD_NUM == 0:
		FOLD_NUM = None
	TRAIN_TEST_RATIO = args.TRAIN_TEST_RATIO
	POS_NEG_RATIO = args.POS_NEG_RATIO
	NUM_EP = args.NUM_EP
	NUM_FRC = args.NUM_FRC
	print("  * K_D_TRAIN: ", K_D_TRAIN)
	print("  * FOLD_NUM: ", FOLD_NUM)
	print("  * TRAIN_TEST_RATIO: ", TRAIN_TEST_RATIO)
	print("  * POS_NEG_RATIO: ", POS_NEG_RATIO)
	print("  * NUM_EP: ", NUM_EP)
	print("  * NUM_FRC: ", NUM_FRC)
	##########################################
	## End of "PASE INPUT"
	##########################################

	##########################################
	## Create feature combination
	##########################################
 	if args.feature_selection == "00000000":
		print "Select at least one feature"
		sys.exit()
	elif args.feature_selection == "000000000":
		this_scores = utils.get_fs_comb(args.feature_selection)
		print "\t".join([fs.name for fs in this_scores])
 	elif args.feature_selection == "000000001":
		this_scores = utils.get_fs_comb(args.feature_selection)
		print "\t".join([fs.name for fs in this_scores])
	else:
		this_scores = utils.get_fs_comb(args.feature_selection)
		print "\t".join([fs.name for fs in this_scores])
	##########################################
	## End of feature combination
	##########################################

	##########################################
	## "Load elution data"
	##########################################
	foundprots, elution_datas = utils.load_data(args.input_dir, this_scores, fc=args.frac_count, mfc=args.elution_max_count)
	##########################################
	## End of "Load elution data"
	##########################################

	##########################################
	## "Generate reference data set"
	##########################################
	gs = ""
	if ((args.taxid != "" and  args.ppi != "") or (args.cluster != "" and  args.ppi != "" )):
		print "Refernce from cluster and PPI are nor compatiple. Please supply ppi or complex reference, not both!"
		sys.exit()

	if args.taxid == "" and  args.ppi == "" and args.cluster == "":
		print "Please supply a reference by setting taxid, cluster, or ppi tag"
		sys.exit()

	gs_clusters = []
	if (args.taxid != "" and args.cluster == "" and args.ppi == ""):
		print "Loading clusters from GO, CORUM, and Intact"
		gs_clusters.extend(utils.get_reference_from_net(args.taxid))

	if args.cluster != "":
		print "Loading complexes from file"
		print("&& args.cluster: ", args.cluster)
		print("&& args.mode: ", args.mode)
		if args.mode == "FA":
			gs_clusters.append(GS.FileClusters(args.cluster, "all"))
		else:
			gs_clusters.append(GS.FileClusters(args.cluster, foundprots))


	if args.ppi != "":
		print "Reading PPI file from %s" % args.reference
		gs = Goldstandard_from_PPI_File(args.ppi, foundprots)

	print gs_clusters
	if 	len(gs_clusters)>0:
		gs = utils.create_goldstandard(gs_clusters, args.taxid, foundprots, pos_neg_ratio = POS_NEG_RATIO)


	output_dir = args.output_dir + os.sep + args.output_prefix

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	refFH = open(output_dir + ".ref_complexes.txt", "w")
	for comp in gs.complexes.complexes:
			print >> refFH, "%s\t%s" % (",".join(comp), ",".join(gs.complexes.complexes[comp]))
	refFH.close()
	##########################################
	## End of "Generate reference data set"
	##########################################


	##########################################
	## "Creating classifier"
	##########################################
	classifier_select = args.classifier
	clf = CS.CLF_Wrapper(args.num_cores, classifier_select, num_ep=NUM_EP, num_frc=NUM_FRC, pos_neg_ratio = POS_NEG_RATIO)
	##########################################
	## End of "Creating classifier"
	##########################################


	##########################################
	## Calculating coelutionScores
	##########################################
	scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=args.num_cores, cutoff= args.co_elution_cutoff)
	if args.precalcualted_score_file == "NONE":
		scoreCalc.calculate_coelutionDatas(gs)
	else:
 		scoreCalc.readTable(args.precalcualted_score_file, gs)

	print("scoreCalc.scores.shape: ", scoreCalc.scores.shape)
	##########################################
	## End of "Calculating coelutionScores"
	##########################################


	##########################################
	## "Balancing Positive & Negative" PPIs
	##########################################
	functionalData = ""
	gs.positive = set(gs.positive & set(scoreCalc.ppiToIndex.keys()))
	gs.negative = set(gs.negative & set(scoreCalc.ppiToIndex.keys()))
	gs.all_positive = gs.positive
	gs.all_negative = gs.negative
	gs.rebalance()

	print("len(gs.positive): ", len(gs.positive))
	print("len(gs.negative): ", len(gs.negative))
	###############################################
	## End of "Balancing Positive & Negative" PPIs
	###############################################

	if args.mode != "EXP":
		print "Loading functional data"
		functionalData = utils.get_FA_data(args.fun_anno_source, args.taxid, args.fun_anno_file)
		print "Dimension of fun anno " + str(functionalData.scores.shape)


	##########################################
	## This part is "EVALUATION"
	##########################################
	print "Start benchmarking"

	if args.mode == "EXP":
		utils.cv_bench_clf(scoreCalc, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO, num_ep = NUM_EP, num_frc = NUM_FRC)

	if args.mode == "COMB":
		tmp_sc = copy.deepcopy(scoreCalc)
		tmp_sc.add_fun_anno(functionalData)
		utils.cv_bench_clf(tmp_sc, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO, num_ep = NUM_EP, num_frc = NUM_FRC)

	if args.mode == "FA":
		utils.cv_bench_clf(functionalData, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO, num_ep = NUM_EP, num_frc = NUM_FRC)
	##########################################
	## End of "EVALUATION"
	##########################################


	# scoreCalc_pkl_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/output/ppi/clf_intensity_H.pkl"


	# scoreCalc_pkl_path = "/Users/chaokuan-hao/Documents/BIO_IT_Station/output/ppi/clf_RAW_1_CNN.pkl"
	# with open(scoreCalc_pkl_path, 'wb') as f:  # Python 3: open(..., 'wb')
	#     pickle.dump(scoreCalc, f)


	# PPI evaluation
	# print utils.cv_bench_clf(scoreCalc, clf, gs, args.output_dir, verbose=False, format="pdf", folds=5)
	#print "I am here"

	##############################################
	## This part is "MODEL TRAINING & PREDICTION"
	##############################################
	network = utils.make_predictions(scoreCalc, args.mode, clf, gs, output_dir, fun_anno=functionalData, verbose = False, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)

	# Predict protein interaction
	outFH = open("%s.pred.txt" % (output_dir), "w")

	final_network = []
	for PPI in network:
		items = PPI.split("\t")
		if float(items[2]) >= args.classifier_cutoff:
		# if float(items[2]) >= 0.0:
			final_network.append(PPI)

	print >> outFH, "\n".join(final_network)
	outFH.close()
	#############################################
	## End of "MODEL TRAINING & PREDICTION"
	#############################################

	# Predicting clusters
	utils.predict_clusters("%s.pred.txt" % (output_dir), "%s.clust.txt" % (output_dir))


	# Evaluating predicted clusters
	pred_clusters = GS.Clusters(False)
	pred_clusters.read_file("%s.clust.txt" % (output_dir))
	overlapped_complexes_with_reference = gs.get_complexes().get_overlapped_complexes_set(pred_clusters)
	print "# of complexes in reference dataset: " + str(len(overlapped_complexes_with_reference))
	#clust_scores, header = utils.clustering_evaluation(gs.complexes, pred_clusters, "", False)
	clust_scores, header, composite_score = utils.clustering_evaluation(gs.complexes, pred_clusters, "", False)
	outFH = open("%s.eval.txt" % (output_dir), "w")
	header = header.split("\t")
	clust_scores = clust_scores.split("\t")
	for i, head in enumerate(header):
		print "%s\t%s" % (head, clust_scores[i])
		print >> outFH, "%s\t%s" % (head, clust_scores[i])
	outFH.close()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

#11000100 (MI, Bayes, PCC+N)
