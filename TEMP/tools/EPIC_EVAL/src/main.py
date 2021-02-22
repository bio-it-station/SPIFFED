from __future__ import division

import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
import sys
import copy
import os
import argparse
import random
import math

import warnings
warnings.filterwarnings('ignore')

def Goldstandard_from_cluster_File(gsF, foundprots=""):
	clusters = GS.Clusters(need_to_be_mapped=False)
	clusters.read_file(gsF)
	if foundprots != "": clusters.remove_proteins(foundprots)
	gs = GS.Goldstandard_from_Complexes("All")
	gs.complexes = clusters
	gs.make_pos_neg_ppis()
	return gs


def Goldstandard_from_PPI_File(gsF, foundprots=""):
	out = GS.Goldstandard_from_Complexes("gs")
	gsFH = open(gsF)
	for line in gsFH:
		line = line.rstrip()
		ida, idb, class_label = line.split("\t")[0:3]
		if foundprots !="" and (ida not in foundprots or idb not in foundprots): continue
		edge = "\t".join(sorted([ida, idb]))
		if class_label == "positive":
			out.positive.add(edge)
		else:
			out.negative.add(edge)
	gsFH.close()
	return out

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
	# n_estimators number:
	parser.add_argument("--N_ESTIMATORS", type = int, help="n_estimators in random forest default = 100", default=100)
	# max_depth number:
	parser.add_argument("--MAX_DEPTH", type = int, help="max_depth in random forest default = None", default=-1)
	# k-fold or direct training:
	parser.add_argument("--K_D_TRAIN", type = str, help="Values: d, k", default="d")
	# fold_num
	parser.add_argument("--FOLD_NUM", type = int, help="default = 5", default=5)
	# train_test_ratio
	parser.add_argument("--TRAIN_TEST_RATIO", type = float, help="default = 0.3", default=0.3)
	# pos_neg_ratio
	parser.add_argument("--POS_NEG_RATIO", type = int, help="default = 1", default=1)
	# eval_ratio
	parser.add_argument("--EVAL_RATIO", type = float, help="default = 0.3", default=0.3)
	############################
	## New added argument END ##
	############################

	args = parser.parse_args()

	args.mode = args.mode.upper()
	args.fun_anno_source = args.fun_anno_source.upper()

	#Create feature combination
 	if args.feature_selection == "00000000":
		print "Select at least one feature"
		sys.exit()

	this_scores = utils.get_fs_comb(args.feature_selection)
	print "\t".join([fs.name for fs in this_scores])

	##########################################
	## This part is "PASE INPUT"
	##########################################
	# N_ESTIMATORS = input("Enter n_estimators:")
	# MAX_DEPTH = input("Enter max_depth:")
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
	# POS_NEG_RATIO = input("Enter negative/positive ratio:")

	N_ESTIMATORS = args.N_ESTIMATORS
	MAX_DEPTH = args.MAX_DEPTH
	if MAX_DEPTH == -1:
		MAX_DEPTH = None
	K_D_TRAIN = args.K_D_TRAIN
	FOLD_NUM = args.FOLD_NUM
	if FOLD_NUM == 0:
		FOLD_NUM = None
	TRAIN_TEST_RATIO = args.TRAIN_TEST_RATIO
	POS_NEG_RATIO = args.POS_NEG_RATIO
	EVAL_RATIO = args.EVAL_RATIO
	print("  * N_ESTIMATORS: ", N_ESTIMATORS)
	print("  * MAX_DEPTH: ", MAX_DEPTH)
	print("  * K_D_TRAIN: ", K_D_TRAIN)
	print("  * FOLD_NUM: ", FOLD_NUM)
	print("  * TRAIN_TEST_RATIO: ", TRAIN_TEST_RATIO)
	print("  * POS_NEG_RATIO: ", POS_NEG_RATIO)
	print("  * EVAL_RATIO: ", EVAL_RATIO)
	##########################################
	## End of "PASE INPUT"
	##########################################

	##########################################
	## "Initialize CLF"
	##########################################
	use_rf = args.classifier == "RF"
	clf = CS.CLF_Wrapper(args.num_cores, use_rf, pos_neg_ratio=POS_NEG_RATIO, n_estimators=N_ESTIMATORS, max_depth=MAX_DEPTH)
	##########################################
	## End of "Initialize CLF"
	##########################################


	##########################################
	## "Load elution data"
	##########################################
 	foundprots, elution_datas = utils.load_data(args.input_dir, this_scores, fc=args.frac_count, mfc=args.elution_max_count)

	print("len(foundprots):", len(foundprots))
	print("len(elution_datas): ", len(elution_datas))
	prots = list(foundprots)
	##########################################
	## End of "Load elution data"
	##########################################

	##########################################
	## This part is "ALL PROTEIN PARTITION"
	##########################################
	prots_len = len(prots)
	split_ratio = EVAL_RATIO
	prots_rand = random.sample(prots, prots_len)
	split_index = int(math.ceil(prots_len*split_ratio))
	prots_eval = prots_rand[0:split_index]
	prots_train = prots_rand[split_index:prots_len]
	foundprots_eval = set(prots_eval)
	foundprots_train = set(prots_train)
	print("len(foundprots_train): ", len(foundprots_train))
	print("len(foundprots_eval): ", len(foundprots_eval))
	##########################################
	## End of "ALL PROTEIN PARTITION"
	##########################################


	##########################################
	## This part is "GENERATE REFERENCE DATA"
	##########################################
	gs = ""
	gs_eval = ""
	if ((args.taxid != "" and  args.ppi != "") or (args.cluster != "" and  args.ppi != "" )):
		print "Refernce from cluster and PPI are nor compatiple. Please supply ppi or complex reference, not both!"
		sys.exit()

	if args.taxid == "" and  args.ppi == "" and args.cluster == "":
		print "Please supply a reference by setting taxid, cluster, or ppi tag"
		sys.exit()

	gs_clusters = []
	gs_eval_clusters = []
	if (args.taxid != "" and args.cluster == "" and args.ppi == ""):
		print "Loading clusters from GO, CORUM, and Intact"
		gs_clusters.extend(utils.get_reference_from_net(args.taxid))
		gs_eval_clusters.extend(utils.get_reference_from_net(args.taxid))

	if args.cluster != "":
		print "Loading complexes from file"
		print("&& args.cluster: ", args.cluster)
		print("&& args.mode: ", args.mode)
		if args.mode == "FA":
			gs_clusters.append(GS.FileClusters(args.cluster, "all"))
		else:
			# gs_clusters.append(GS.FileClusters(args.cluster, foundprots))
			# Goes here (change from foundprots to foundprots_final)
			gs_clusters.append(GS.FileClusters(args.cluster, foundprots_train))
			gs_eval_clusters.append(GS.FileClusters(args.cluster, foundprots_eval))

	if args.ppi != "":
		print "Reading PPI file from %s" % args.reference
		# gs = Goldstandard_from_PPI_File(args.ppi, foundprots)
		gs = Goldstandard_from_PPI_File(args.ppi, foundprots_train)
		gs_eval = Goldstandard_from_PPI_File(args.ppi, foundprots_eval)


	print("gs_clusters: ", gs_clusters)
	print("gs_eval_clusters: ", gs_eval_clusters)
	if 	len(gs_clusters)>0:
		# Goes here (change from foundprots to foundprots_final)
		print("Training utils.create_goldstandard: ")
		gs = utils.create_goldstandard(gs_clusters, args.taxid, foundprots, pos_neg_ratio = POS_NEG_RATIO)

	if 	len(gs_eval_clusters)>0:
		# Goes here (change from foundprots to foundprots_eval_final)
		print("Testing utils.create_goldstandard: ")
		gs_eval = utils.create_goldstandard(gs_eval_clusters, args.taxid, foundprots_eval, pos_neg_ratio = POS_NEG_RATIO)
	##########################################
	## End of "GENERATE REFERENCE DATA"
	##########################################


	##########################################
	## This part is "WRITING GS COMPLEXES"
	##########################################
	output_dir = args.output_dir + os.sep + args.output_prefix
	# print("output_dir: ", output_dir)

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	refFH = open(output_dir + ".ref_complexes.txt", "w")
	for comp in gs.complexes.complexes:
			print >> refFH, "%s\t%s" % (",".join(comp), ",".join(gs.complexes.complexes[comp]))
	refFH.close()

	refFH_eval = open(output_dir + ".ref_eval_complexes.txt", "w")
	for comp in gs_eval.complexes.complexes:
		print >> refFH_eval, "%s\t%s" % (",".join(comp), ",".join(gs_eval.complexes.complexes[comp]))
	refFH_eval.close()
	##########################################
	## End of "WRITING GS COMPLEXES"
	##########################################


	##########################################
	## This part is "CALCULATE ELUTION SCORES"
	##########################################
	scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=args.num_cores, cutoff= args.co_elution_cutoff)
	scoreCalc_eval = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".eval_scores.txt", num_cores=args.num_cores, cutoff= args.co_elution_cutoff)

	if args.precalcualted_score_file == "NONE":
		scoreCalc.calculate_coelutionDatas(gs, list(foundprots_train))
		scoreCalc_eval.calculate_coelutionDatas(gs_eval, list(foundprots_eval))
	else:
 		scoreCalc.readTable(args.precalcualted_score_file, gs)
	 	scoreCalc_eval.readTable(args.precalcualted_score_file, gs_eval)

	print("scoreCalc.scores.shape: ", scoreCalc.scores.shape)
	print("scoreCalc_eval.scores.shape: ", scoreCalc_eval.scores.shape)
	##########################################
	## End of "CALCULATE ELUTION SCORES"
	##########################################

	##########################################
	## "Balancing Positive & Negative" PPIs
	##########################################
	# Remove +/- PPIs that are removed
	## This is for training
	functionalData = ""
	gs.positive = set(gs.positive & set(scoreCalc.ppiToIndex.keys()))
	gs.negative = set(gs.negative & set(scoreCalc.ppiToIndex.keys()))
	gs.all_positive = gs.positive
	gs.all_negative = gs.negative
	gs.rebalance()

	print("len(gs.positive): ", len(gs.positive))
	print("len(gs.negative): ", len(gs.negative))
	print("len(gs.all_positive): ", len(gs.all_positive))
	print("len(gs.all_negative): ", len(gs.all_negative))

	## This is for testing
	functionalEvalData = ""
	gs_eval.positive = set(gs_eval.positive & set(scoreCalc_eval.ppiToIndex.keys()))
	gs_eval.negative = set(gs_eval.negative & set(scoreCalc_eval.ppiToIndex.keys()))
	gs_eval.all_positive = gs_eval.positive
	gs_eval.all_negative = gs_eval.negative
	gs_eval.rebalance()

	print("len(gs_eval.positive): ", len(gs_eval.positive))
	print("len(gs_eval.negative): ", len(gs_eval.negative))
	print("len(gs_eval.all_positive): ", len(gs_eval.all_positive))
	print("len(gs_eval.all_negative): ", len(gs_eval.all_negative))
	###############################################
	## End of "Balancing Positive & Negative" PPIs
	###############################################

	# Not entry
	if args.mode != "EXP":
		print "Loading functional data"
		functionalData = utils.get_FA_data(args.fun_anno_source, args.taxid, args.fun_anno_file)
		print "Dimension of fun anno " + str(functionalData.scores.shape)


	##########################################
	## This part is "EVALUATION"
	##########################################
	# if K_D_TRAIN == "k":
	# 	pass
	# elif K_D_TRAIN == "d":
	# print "Start benchmarking"
	#
	# if args.mode == "EXP":
	# 	utils.cv_bench_clf(scoreCalc, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)
	#
	# if args.mode == "COMB":
	# 	tmp_sc = copy.deepcopy(scoreCalc)
	# 	tmp_sc.add_fun_anno(functionalData)
	# 	utils.cv_bench_clf(tmp_sc, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)
	#
	# if args.mode == "FA":
	# 	utils.cv_bench_clf(functionalData, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)

	# PPI evaluation
	# print utils.cv_bench_clf(scoreCalc, clf, gs, args.output_dir, verbose=False, format="pdf", k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)
	##########################################
	## End of "EVALUATION"
	##########################################


	##########################################
	## This part is "MODEL TRAINING & PREDICTION"
	##########################################
	# network = utils.make_predictions(scoreCalc, args.mode, clf, gs, output_dir, fun_anno=functionalData, verbose = False, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)
	network, network_eval = utils.make_predictions(scoreCalc, scoreCalc_eval, args.mode, clf, gs, gs_eval, output_dir, fun_anno=functionalData, verbose = False, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)

	# Predict protein interaction
	outFH = open("%s.pred.txt" % (output_dir), "w")


	final_network = []
	for PPI in network:
		items = PPI.split("\t")
		if float(items[2]) >= args.classifier_cutoff:
		# if float(items[2]) >= 0.0:
			final_network.append(PPI)
	print >> outFH, "\n".join(final_network)

	final_network_eval = []
	for PPI in network_eval:
		items = PPI.split("\t")
		if float(items[2]) >= args.classifier_cutoff:
		# if float(items[2]) >= 0.0:
			final_network_eval.append(PPI)
	print >> outFH, "\n".join(final_network_eval)


	outFH.close()
	##########################################
	## End of "MODEL TRAINING & PREDICTION"
	##########################################

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
