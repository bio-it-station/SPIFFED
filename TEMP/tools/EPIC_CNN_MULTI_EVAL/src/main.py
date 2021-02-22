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
import random
import math

import warnings
warnings.filterwarnings('ignore')

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

	args = parser.parse_args()
	args.mode = args.mode.upper()
	args.fun_anno_source = args.fun_anno_source.upper()


	##########################################
	## This part is "PASE INPUT"
	##########################################
	K_D_TRAIN = raw_input("Enter k-fold or direct training (k/d):")
	if K_D_TRAIN == "k":
		print "You are using k-fold training!"
		FOLD_NUM = input("Enter fold number:")
		print("Fold number is: " + str(FOLD_NUM))
		TRAIN_TEST_RATIO = None
	elif K_D_TRAIN == "d":
		print "You are using direct training!"
		TRAIN_TEST_RATIO = input("Enter test/all_data ratio:")
		print("Your TRAIN_TEST_RATIO ratio is: " + str(TRAIN_TEST_RATIO))
		FOLD_NUM = None
	else:
		print "Please input 'k' or 'd'."
		sys.exit()
	NUM_EP = 2
	NUM_FRC = input("Enter number of fractions:")
	POS_NEG_RATIO = input("Enter negative/positive ratio:")
	##########################################
	## End of "PASE INPUT"
	##########################################


	##########################################
	## This part is "FEATURE COMB CREATION"
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
	## End of "FEATURE COMB CREATION"
	##########################################


	##########################################
	## This part is "GETTING ALL PROTEIN"
	##########################################
	input_dirs = []
	for each in os.listdir(args.input_dir):
		for file in os.listdir(os.path.join(args.input_dir, each)):
			input_dirs.append(os.path.join(args.input_dir, each, file))

	print "input_dirs: ", input_dirs

	foundprots_ens, elution_datas_ens = utils.load_data(input_dirs, this_scores, fc=args.frac_count, mfc=args.elution_max_count)
	# print("foundprots_ens: ", foundprots_ens)
	print("len(foundprots_ens):", len(foundprots_ens))
	# print("elution_datas_ens: ", elution_datas_ens)
	print("len(elution_datas_ens): ", len(elution_datas_ens))

	prots_ens = list(foundprots_ens)
	##########################################
	## End of "GETTING ALL PROTEIN"
	##########################################


	##########################################
	## This part is "ALL PROTEIN PARTITION"
	##########################################
	prots_ens_len = len(prots_ens)
	split_ratio = TRAIN_TEST_RATIO
	prots_ens_rand = random.sample(prots_ens, prots_ens_len)
	split_index = int(math.ceil(prots_ens_len*split_ratio))
	prots_ens_eval = prots_ens_rand[0:split_index]
	prots_ens_train = prots_ens_rand[split_index:prots_ens_len]

	# print("prots_ens_train: ", prots_ens_train)
	print("len(prots_ens_train): ", len(prots_ens_train))
	# print("prots_ens_eval: ", prots_ens_eval)
	print("len(prots_ens_eval): ", len(prots_ens_eval))
	##########################################
	## End of "ALL PROTEIN PARTITION"
	##########################################














	networks = []
	networks_eval = []
	PPI_voting_dic = {}
	PPI_voting_eval_dic = {}
	print("args.input_dir: ", args.input_dir)
	for each in os.listdir(args.input_dir):
		input_dir = os.path.join(args.input_dir, each)
		print("each: ", each)
		# Load elution data
		##########################################
		## This part is "LOAD ELUTION DATA"
		##########################################
		# I need to find intersection between foundprots and prots_ens_train
	 	foundprots, elution_datas = utils.load_data(input_dir, this_scores, fc=args.frac_count, mfc=args.elution_max_count)

		## Here is the training data for ensemble model
		foundprots_train = set(prots_ens_train)
		foundprots_final = foundprots_train.intersection(foundprots)

		## Here is the testing data for ensemble model
		foundprots_eval = set(prots_ens_eval)
		foundprots_eval_final = foundprots_eval.intersection(foundprots)

		print("len(foundprots): ", len(foundprots))
		print("len(foundprots_train): ", len(foundprots_train))
		print("len(foundprots_final): ", len(foundprots_final))
		print("len(foundprots_eval): ", len(foundprots_eval))
		print("len(foundprots_eval_final): ", len(foundprots_eval_final))
		##########################################
		## End of "LOAD ELUTION DATA"
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
				# Goes here (change from foundprots to foundprots_final)
				gs_clusters.append(GS.FileClusters(args.cluster, foundprots_final))
				gs_eval_clusters.append(GS.FileClusters(args.cluster, foundprots_eval_final))

		if args.ppi != "":
			print "Reading PPI file from %s" % args.reference
			gs = Goldstandard_from_PPI_File(args.ppi, foundprots_final)
			gs_eval = Goldstandard_from_PPI_File(args.ppi, foundprots_eval_final)

		print("gs_clusters: ", gs_clusters)
		print("gs_eval_clusters: ", gs_eval_clusters)
		if 	len(gs_clusters)>0:
			# Goes here (change from foundprots to foundprots_final)
			print("Training utils.create_goldstandard: ")
			gs = utils.create_goldstandard(gs_clusters, args.taxid, foundprots_final, pos_neg_ratio = POS_NEG_RATIO)

		if 	len(gs_eval_clusters)>0:
			# Goes here (change from foundprots to foundprots_eval_final)
			print("Testing utils.create_goldstandard: ")
			gs_eval = utils.create_goldstandard(gs_eval_clusters, args.taxid, foundprots_eval_final, pos_neg_ratio = POS_NEG_RATIO)
		##########################################
		## End of "GENERATE REFERENCE DATA"
		##########################################


		##########################################
		## This part is "WRITING GS COMPLEXES"
		##########################################
		output_dir = args.output_dir + os.sep + each + os.sep +  args.output_prefix

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
		## This part is "CLASSIFIER"
		##########################################
		classifier_select = args.classifier
		clf = CS.CLF_Wrapper(args.num_cores, classifier_select, num_ep=NUM_EP, num_frc=NUM_FRC)
		##########################################
		## End of "CLASSIFIER"
		##########################################


		##########################################
		## This part is "CALCULATE ELUTION SCORES"
		##########################################
		scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=args.num_cores, cutoff= args.co_elution_cutoff)
		scoreCalc_eval = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".eval_scores.txt", num_cores=args.num_cores, cutoff= args.co_elution_cutoff)

		if args.precalcualted_score_file == "NONE":
			scoreCalc.calculate_coelutionDatas(gs, list(foundprots_final))
			scoreCalc_eval.calculate_coelutionDatas(gs_eval, list(foundprots_eval_final))
		else:
	 		scoreCalc.readTable(args.precalcualted_score_file, gs)
	 		scoreCalc_eval.readTable(args.precalcualted_score_file, gs_eval)
		print("scoreCalc.scores.shape: ", scoreCalc.scores.shape)
		print("scoreCalc_eval.scores.shape: ", scoreCalc_eval.scores.shape)
		##########################################
		## End of "CALCULATE ELUTION SCORES"
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

		## This is for testing
		functionalEvalData = ""
		gs_eval.positive = set(gs_eval.positive & set(scoreCalc_eval.ppiToIndex.keys()))
		gs_eval.negative = set(gs_eval.negative & set(scoreCalc_eval.ppiToIndex.keys()))
		gs_eval.all_positive = gs_eval.positive
		gs_eval.all_negative = gs_eval.negative
		gs_eval.rebalance()

		print("len(gs_eval.positive): ", len(gs_eval.positive))
		print("len(gs_eval.negative): ", len(gs_eval.negative))


		# Not entry
		if args.mode != "EXP":
			print "Loading functional data"
			functionalData = utils.get_FA_data(args.fun_anno_source, args.taxid, args.fun_anno_file)
			print "Dimension of fun anno " + str(functionalData.scores.shape)


		# ##########################################
		# ## This part is "EVALUATION"
		# ##########################################
		# print "Start benchmarking"
		#
		# if args.mode == "EXP":
		# 	utils.cv_bench_clf(scoreCalc, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO, num_ep = NUM_EP, num_frc = NUM_FRC)
		#
		# if args.mode == "COMB":
		# 	tmp_sc = copy.deepcopy(scoreCalc)
		# 	tmp_sc.add_fun_anno(functionalData)
		# 	utils.cv_bench_clf(tmp_sc, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO, num_ep = NUM_EP, num_frc = NUM_FRC)
		#
		# if args.mode == "FA":
		# 	utils.cv_bench_clf(functionalData, clf, gs, output_dir, format="pdf", verbose=True, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO, num_ep = NUM_EP, num_frc = NUM_FRC)
		# ##########################################
		# ## End of "EVALUATION"
		# ##########################################


		##########################################
		## This part is "MODEL TRAINING & PREDICTION"
		##########################################
		network, network_eval = utils.make_predictions(scoreCalc, scoreCalc_eval, args.mode, clf, gs, gs_eval, output_dir, fun_anno=functionalData, verbose = False, k_d_training = K_D_TRAIN, folds = FOLD_NUM, train_test_ratio = TRAIN_TEST_RATIO)
		## Adding network into the networks list
		networks.append(network)
		networks_eval.append(network_eval)

		# Predict protein interaction (Training)
		outFH = open("%s.pred.txt" % (output_dir), "w")
		final_network = []
		for PPI in network:
			items = PPI.split("\t")
			if float(items[2]) >= args.classifier_cutoff:
				####################################
				## This part is "ENSEMBLE PPI CREAT"
				####################################
				key = items[0]+'\t'+items[1]
				if key in PPI_voting_dic:
					PPI_voting_dic[key].append(float(items[2]))
				else:
					PPI_voting_dic[key] = [float(items[2])]
				####################################
				## End of "ENSEMBLE PPI CREAT"
				####################################
				final_network.append(PPI)

		print >> outFH, "\n".join(final_network)
		outFH.close()


		# Predict protein interaction (Testing)
		outFH_eval = open("%s.eval.pred.txt" % (output_dir), "w")
		final_network_eval = []
		for PPI in network_eval:
			items = PPI.split("\t")
			if float(items[2]) >= args.classifier_cutoff:
				####################################
				## This part is "ENSEMBLE PPI CREAT"
				####################################
				key = items[0]+'\t'+items[1]
				if key in PPI_voting_eval_dic:
					PPI_voting_eval_dic[key].append(float(items[2]))
				else:
					PPI_voting_eval_dic[key] = [float(items[2])]
				####################################
				## End of "ENSEMBLE PPI CREAT"
				####################################
				final_network_eval.append(PPI)

		print >> outFH_eval, "\n".join(final_network_eval)
		outFH_eval.close()
		##########################################
		## End of "MODEL TRAINING & PREDICTION"
		##########################################


	##########################################
	## This part is "ENSEMBLE VOTING"
	##########################################
	strict_majority = (len(networks)+2-1) // 2

	# This is the normal voting
	output_dir_ens = args.output_dir + os.sep + "vote" + os.sep + args.output_prefix

	outdir_str_maj = output_dir_ens + "_str_maj.pred.txt"
	outdir_maj = output_dir_ens + "_maj.pred.txt"
	outdir_one = output_dir_ens + "_1.pred.txt"

	final_network_str_maj = []
	final_network_maj = []
	final_network_one = []

	with open(outdir_str_maj, "w") as f1, open(outdir_maj, "w") as f2, open(outdir_one, "w") as f3:
		for PPI, edges in PPI_voting_dic.iteritems():
			edges = np.array(edges)
			avg_pos_edge = np.average(edges[edges >= 0.5])
			ppi_edges = PPI + '\t' + str(avg_pos_edge)

			# This is the condition for strict majority vote
			pos_count = np.sum(edges >= 0.5)
			if pos_count >= strict_majority:
				final_network_str_maj.append(ppi_edges)

			# This is the condition for majority vote
			majority = (len(edges)+2-1) // 2
			if pos_count >= majority:
				final_network_maj.append(ppi_edges)

			# This is the condition for one-agreement vote
			any_pos = any(edges >= 0.5)
			if any_pos:
				final_network_one.append(ppi_edges)

		print >> f1, "\n".join(final_network_str_maj)
		print >> f2, "\n".join(final_network_maj)
		print >> f3, "\n".join(final_network_one)


	# This is the voting evaluation
	output_dir_eval_ens = args.output_dir + os.sep + "vote_eval" + os.sep + args.output_prefix

	outdir_str_maj_eval = output_dir_eval_ens + "_str_maj.pred.txt"
	outdir_maj_eval = output_dir_eval_ens + "_maj.pred.txt"
	outdir_one_eval = output_dir_eval_ens + "_1.pred.txt"

	final_network_str_maj_eval = []
	final_network_maj_eval = []
	final_network_one_eval = []

	with open(outdir_str_maj_eval, "w") as f1, open(outdir_maj_eval, "w") as f2, open(outdir_one_eval, "w") as f3:
		for PPI, edges in PPI_voting_eval_dic.iteritems():
			edges = np.array(edges)
			avg_pos_edge = np.average(edges[edges >= 0.5])
			ppi_edges = PPI + '\t' + str(avg_pos_edge)

			# This is the condition for strict majority vote
			pos_count = np.sum(edges >= 0.5)
			if pos_count >= strict_majority:
				final_network_str_maj_eval.append(ppi_edges)

			# This is the condition for majority vote
			majority = (len(edges)+2-1) // 2
			if pos_count >= majority:
				final_network_maj_eval.append(ppi_edges)

			# This is the condition for one-agreement vote
			any_pos = any(edges >= 0.5)
			if any_pos:
				final_network_one_eval.append(ppi_edges)

		print >> f1, "\n".join(final_network_str_maj_eval)
		print >> f2, "\n".join(final_network_maj_eval)
		print >> f3, "\n".join(final_network_one)
	##########################################
	## End of "ENSEMBLE VOTING"
	##########################################


	##########################################
	## This part is "Reference Complex Creat"
	##########################################
	gs_clusters_ens = []
	gs_clusters_ens.append(GS.FileClusters(args.cluster, foundprots_ens))
	if len(gs_clusters_ens)>0:
		gs_ens = utils.create_goldstandard(gs_clusters_ens, args.taxid, foundprots_ens, pos_neg_ratio = POS_NEG_RATIO)

	refFH = open(output_dir_ens + ".ref_complexes.txt", "w")
	for comp in gs_ens.complexes.complexes:
		print >> refFH, "%s\t%s" % (",".join(comp), ",".join(gs_ens.complexes.complexes[comp]))
	refFH.close()

	gs_ens.all_positive = gs_ens.positive
	gs_ens.all_negative = gs_ens.negative
	gs_ens.rebalance()



	##########################################
	## End of "Reference Complex Creat"
	##########################################

	for suffix in ["_str_maj", "_maj", "_1"]:
		# Predicting clusters
		utils.predict_clusters(output_dir_ens+suffix+".pred.txt", output_dir_ens+suffix+".clust.txt")

		# Evaluating predicted clusters
		pred_clusters = GS.Clusters(False)
		pred_clusters.read_file(output_dir_ens+suffix+".clust.txt")
		overlapped_complexes_with_reference = gs_ens.get_complexes().get_overlapped_complexes_set(pred_clusters)
		print "This is :" + suffix
		print "# of complexes in reference dataset: " + str(len(overlapped_complexes_with_reference))
		clust_scores, header, composite_score = utils.clustering_evaluation(gs_ens.complexes, pred_clusters, "", False)

		with open(output_dir_ens+suffix+".eval.txt", "w") as f:
		# outFH = open(, "w")
			header = header.split("\t")
			clust_scores = clust_scores.split("\t")
			for i, head in enumerate(header):
				print "%s\t%s" % (head, clust_scores[i])
				print >> f, "%s\t%s" % (head, clust_scores[i])
		# outFH.close()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

#11000100 (MI, Bayes, PCC+N)
