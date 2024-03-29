from __future__ import division

import CalculateCoElutionScores as CS
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import GoldStandard as GS
import copy
import json

from keras.callbacks import EarlyStopping
from sklearn.model_selection import train_test_split
from sklearn import metrics
import scikitplot as skplt
import numpy as np


def get_fs_comb(comb_string):
	#Create feature combination
	# scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	# this_scores = []
	# for i, feature_selection in enumerate(comb_string):
	# 	if feature_selection == "1": this_scores.append(scores[i])
	# return this_scores
	#Create feature combination
	if comb_string == "000000000":
		this_scores = [CS.LocalPcc()]
	else:
		scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex(), CS.Raw_eps()]
		this_scores = []
		for i, feature_selection in enumerate(comb_string):
			if feature_selection == "1": this_scores.append(scores[i])
		print("this_scores: ", this_scores)
	return this_scores

# a new function added by Lucas HU for benchmark on PPIs level
# this is used for two-levels optimization steps
# a trial version
def bench_by_PPI_clf(num_folds, scoreCalc, train_gold_complexes, clf):
	_, data_train, targets_train = scoreCalc.toSklearnData(train_gold_complexes)

	# define the correlation score matrix for positive PPIs and negative PPIs.
	positive_data = np.zeros((sum(targets_train), np.shape(data_train)[1]))
	negative_data = np.zeros((len(targets_train) - sum(targets_train), np.shape(data_train)[1]))

	index = 0
	positive_index = 0
	negative_index = 0

	for label in targets_train:
		if label == 1:
			positive_data[positive_index, :] = data_train[index, :]
			positive_index = positive_index + 1
		if label == 0:
			negative_data[negative_index, :] = data_train[index, :]
			negative_index = negative_index + 1
		index = index + 1

	fold_size_positive = int(np.shape(positive_data)[0] / num_folds)
	fold_size_negative = int(np.shape(negative_data)[0] / num_folds)

	# set the initial values for the three metrics
	fmeasure_sum = 0
	auc_pr_sum = 0
	auc_roc_sum = 0
	precision_sum = 0
	recall_sum = 0


	# do n_fold_cross_validation and reported the avaergae value of all measurement metrics
	for i in range(num_folds):
		eval_positive = positive_data[fold_size_positive * i : fold_size_positive * (i + 1),:]
		index_rows_for_eval_positive = list(range(fold_size_positive * i , fold_size_positive * (i + 1)))
		train_positive = np.delete(positive_data, index_rows_for_eval_positive, 0)

		eval_negative = negative_data[fold_size_negative * i : fold_size_negative * (i + 1),:]
		index_rows_for_eval_negative = list(range(fold_size_negative * i, fold_size_negative * (i + 1)))
		train_negative = np.delete(negative_data, index_rows_for_eval_negative, 0)

		eval_data = np.concatenate((eval_positive, eval_negative), axis=0)
		train_data = np.concatenate((train_positive, train_negative), axis=0)

		eval_positive_labels = np.array([1] * np.shape(eval_positive)[0])
		eval_negative_labels = np.array([0] * np.shape(eval_negative)[0])

		train_positive_labels = np.array([1] * np.shape(train_positive)[0])
		train_negative_labels = np.array([0] * np.shape(train_negative)[0])

		eval_labels = np.concatenate([eval_positive_labels, eval_negative_labels])
		train_labels = np.concatenate([train_positive_labels, train_negative_labels])

		#train the classifier
		clf.fit(train_data, train_labels)

		#evaluate the classifier
		precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.eval(eval_data, eval_labels)
		fmeasure_sum = fmeasure_sum + fmeasure
		auc_pr_sum = auc_pr_sum + auc_pr
		auc_roc_sum = auc_roc_sum + auc_roc
		precision_sum = precision_sum + precision
		recall_sum = recall_sum + recall

		recall_vals, precision_vals, threshold = curve_pr
		threshold = np.append(threshold, 1)

	fmeasure_average = fmeasure_sum / num_folds
	auc_pr_average = auc_pr_sum / num_folds
	auc_roc_average = auc_roc_sum / num_folds
	precision_avergae = precision/ num_folds
	recall_average = recall_sum/ num_folds

	avergae_list = [precision_avergae, recall_average, fmeasure_average, auc_pr_average, auc_roc_average]

	return avergae_list



#@ Kuan-Hao Chao
# Add 'CNN' and 'Label Spreading' methods
def cv_bench_clf(scoreCalc, clf, gs, outDir, verbose=False, learning_selection = 'sl', format="pdf", folds=None, train_test_ratio=None, num_ep = 2, num_frc = None):
# def cv_bench_clf(scoreCalc, clf, gs, outDir, verbose=False, format="pdf", folds = 5):
	_, data, targets, unsure_data, unsure_targets = scoreCalc.toSklearnData(gs)

	print("data: ", data.shape)
	print("data.shape: ", data.shape)
	print("targets.shape: ", targets.shape)

	#######################################
	## "Data reshape"
	#######################################
	if clf.classifier_select == "CNN":
		data = np.array(data)
		num_samples, num_scores, num_frc = data.shape
		print("data: ", data.shape)
		data = data.reshape((num_samples, num_scores, num_frc, 1))
		print("data.shape: ", data.shape)
		print(data[0])
		print("targets.shape: ", targets.shape)
		print(targets)

		unsure_data = np.array(unsure_data)
		num_samples_uns, num_scores_uns, num_frc_uns = unsure_data.shape
		print("unsure_data: ", unsure_data.shape)
		unsure_data = unsure_data.reshape((num_samples_uns, num_scores_uns, num_frc_uns, 1))
		print("unsure_data.shape: ", unsure_data.shape)
		num_ppi_uns, num_ep_uns, num_frc_uns, tmp = unsure_data.shape
	elif clf.classifier_select == "LS":
		data = np.array(data)
		num_samples, num_scores, num_frc = data.shape
		print("data: ", data.shape)
		data = data.reshape((num_samples, num_scores*num_frc))
		print("data.shape: ", data.shape)
		print(data[0])
		print("targets.shape: ", targets.shape)
		print(targets)

		unsure_data = np.array(unsure_data)
		num_samples_uns, num_scores_uns, num_frc_uns = unsure_data.shape
		print("unsure_data: ", unsure_data.shape)
		unsure_data = unsure_data.reshape((num_samples_uns, num_scores_uns*num_frc_uns))
		print("unsure_data.shape: ", unsure_data.shape)
	#######################################
	## End of "Data reshape"
	#######################################


	if learning_selection == 'sl':
		## This can be set to baseline model
		this_targets_test, preds_test, probs_test, precision_test, recall_test, fmeasure_test, auc_pr_test, auc_roc_test, curve_pr_test, curve_roc_test, this_targets_train, preds_train, probs_train, precision_train, recall_train, fmeasure_train, auc_pr_train, auc_roc_train, curve_pr_train, curve_roc_train = clf.cv_eval(data, targets, unsure_data, unsure_targets, outDir, folds, train_test_ratio, num_ep, num_frc)
		outDir_all_pos_neg_testing = outDir + os.sep + "all_pos_neg_testing"
		outDir_all_pos_neg_training = outDir + os.sep + "all_pos_neg_training"

		eval_plotting(outDir_all_pos_neg_testing, this_targets_test, preds_test, probs_test, precision_test, recall_test, fmeasure_test, auc_pr_test, auc_roc_test, curve_pr_test, curve_roc_test)
		eval_plotting(outDir_all_pos_neg_training, this_targets_train, preds_train, probs_train, precision_train, recall_train, fmeasure_train, auc_pr_train, auc_roc_train, curve_pr_train, curve_roc_train)


	if learning_selection == 'ssl':
		## This is the evaluation for semi-supervised learning
		this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.cv_ssl_eval(data, targets, unsure_data, unsure_targets, outDir, folds, train_test_ratio, num_ep, num_frc)
		outDir_all_pos_neg_testing = outDir + os.sep + "ssl_all_pos_neg_testing"
		eval_plotting(outDir_all_pos_neg_testing, this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc)

	rownames = ["Precision", "Recall", "F-Measure", "AUC PR", "AUC ROC"]
	# return rownames, [precision, recall, fmeasure, auc_pr, auc_roc]

# def bench_clf(scoreCalc, train, eval, clf, outDir, verbose=False, format = "pdf"):
# 	_, data_train, targets_train = scoreCalc.toSklearnData(train)
# 	_, data_eval, targets_eval = scoreCalc.toSklearnData(eval)
#
# 	clf.fit(data_train, targets_train)
# 	precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.eval(data_eval, targets_eval)
# 	plotCurves([("", curve_roc)], outDir + ".roc." + format, "False Positive rate", "True Positive Rate")
# 	recall_vals, precision_vals, threshold = curve_pr
# 	plotCurves([("", (precision_vals, recall_vals))], outDir + ".pr." + format, "Recall", "Precision")
#
# 	threshold = np.append(threshold, 1)
# 	plotCurves([("Precision", (precision_vals, threshold)), ("Recall", (recall_vals, threshold))], outDir + ".cutoff." + format, "Cutoff", "Evaluation metric score")
# 	if verbose:
# 		rownames = ["Precision", "Recall", "F-Measure", "AUC PR", "AUC ROC"]
# 		val_scores = [precision, recall, fmeasure, auc_pr, auc_roc]
# 		for i in range(len(rownames)):
# 			print rownames[i]
# 			print val_scores[i]



# a function added by Lucas HU for n_fold corss validation
# a trial verison
def make_predictions_cross_validation(scoreCalc, train, eval, clf):
	_, data_train, targets_train = scoreCalc.toSklearnData(train)
	networkDic = set([])

	eval_names, data_eval, targets_eval = scoreCalc.toSklearnData(eval)
	if len(eval_names) == 0: return networkDic

	print "To pred"
	print data_eval.shape

	tmp_clf = copy.deepcopy(clf)
	tmp_clf.fit(data_train, targets_train)
	probs, predicts = tmp_clf.predict_proba(data_eval), tmp_clf.predict(data_eval)
	for index in range(len(probs)):
		if predicts[index] == 1:
			networkDic.add("%s\t%f" % (eval_names[index], probs[index]))
	return networkDic

# @author: Florian Goebels
# makes precision recall plot for mutliple rp ccurves
# @Param:
#	curves list of tuples with (name, precision, recall) which should be plotted
#	outF Pdf file location for the created plot
def plotCurves(curves, outF, xlab, ylab):
	plt.clf()
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	cols = ['b', 'r', 'c', 'm', 'y', 'k']
	for (name, curve) in curves:
		x, y = curve[0:2]
		if name != "":
			plt.plot(x, y, label=name, color = cols.pop())
		else:
			plt.plot(x, y, color=cols.pop())
	art = []
	if len(curves)>1:
		lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1),  ncol = 5, fontsize=8)
		art.append(lgd)
	plt.savefig(outF, additional_artists=art, bbox_inches="tight")
	plt.close()






# # @author Florian Goebels
# def predictInteractions(scoreCalc, clf, gs, outDir, to_train=True, verbose= True, k_d_training = "k", folds = 5, train_test_ratio = None):
# 	this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.cv_model_creation(data_train, targets_train, k_d_training, folds, train_test_ratio)
#
#
# 	outDir_all_pos_neg_train = outDir + os.sep + "all_pos_neg_training"
# 	eval_plotting(outDir_all_pos_neg_train, this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc)
#
#
#
# 	this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.cv_model_all_pos_neg_eval(data_all, targets_all, k_d_training, folds, train_test_ratio)
#
# 	outDir_all_pos_neg = outDir + os.sep + "all_pos_neg"
# 	eval_plotting(outDir_all_pos_neg, this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc)
#
#
# 	num_features = data_train.shape[1]
#
# 	print("data_train: ", data_train.shape)
# 	print("data_train.shape: ", data_train.shape)
#@ Kuan-Hao Chao
# Add 'CNN' and 'Label Spreading' methods
def predictInteractions(scoreCalc, clf, gs, outDir, to_train=True, verbose= True, folds = 5, train_test_ratio = None, num_ep=2, num_frc=27):

	_, data, targets, unsure_data, unsure_targets = scoreCalc.toSklearnData(gs)
	_all, data_all, targets_all, unsure_data, unsure_targets = scoreCalc.toSklearnDataAll(gs)

	if clf.classifier_select == "RF":
		pass
	#######################################
	## "Data reshape"
	#######################################
	elif clf.classifier_select == "CNN":
		data = np.array(data)
		num_samples, num_scores, num_frc = data.shape
		print("data: ", data.shape)
		data = data.reshape((num_samples, num_scores, num_frc, 1))
		print("data.shape: ", data.shape)
		num_ppi, num_ep, num_frc, tmp = data.shape

		data_all = np.array(data_all)
		num_samples_all, num_scores_all, num_frc_all = data_all.shape
		print("data: ", data_all.shape)
		data_all = data_all.reshape((num_samples_all, num_scores_all, num_frc_all, 1))
		print("data_all.shape: ", data_all.shape)
		num_ppi_all, num_ep_all, num_frc_all, tmp = data_all.shape


		unsure_data = np.array(unsure_data)
		num_samples_uns, num_scores_uns, num_frc_uns = unsure_data.shape
		print("unsure_data: ", unsure_data.shape)
		unsure_data = unsure_data.reshape((num_samples_uns, num_scores_uns, num_frc_uns, 1))
		print("unsure_data.shape: ", unsure_data.shape)
		num_ppi_uns, num_ep_uns, num_frc_uns, tmp = unsure_data.shape
	elif clf.classifier_select == "LS":
		data = np.array(data)
		num_samples, num_scores, num_frc = data.shape
		print("data: ", data.shape)
		data = data.reshape((num_samples, num_scores*num_frc))
		print("data.shape: ", data.shape)
		num_ppi, num_ep_frc = data.shape

		data_all = np.array(data_all)
		num_samples_all, num_scores_all, num_frc_all = data_all.shape
		print("data: ", data_all.shape)
		data_all = data_all.reshape((num_samples_all, num_scores_all*num_frc_all))
		print("data_all.shape: ", data_all.shape)
		num_ppi_all, num_ep_frc_all = data_all.shape

		unsure_data = np.array(unsure_data)
		num_samples_uns, num_scores_uns, num_frc_uns = unsure_data.shape
		print("unsure_data: ", unsure_data.shape)
		unsure_data = unsure_data.reshape((num_samples_uns, num_scores_uns*num_frc_uns))
		print("unsure_data.shape: ", unsure_data.shape)
		num_ppi_uns, num_ep_frc_uns = unsure_data.shape
	#######################################
	## End of "Data reshape"
	#######################################

	this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.cv_model_creation(data, targets, unsure_data, unsure_targets, folds, train_test_ratio, num_ep, num_frc)

	outDir_all_pos_neg_train = outDir + os.sep + "all_pos_neg_training_whole"
	eval_plotting(outDir_all_pos_neg_train, this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc)

	this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.cv_model_all_pos_neg_eval(data_all, targets_all, folds, train_test_ratio)

	outDir_all_pos_neg = outDir + os.sep + "all_pos_neg"
	eval_plotting(outDir_all_pos_neg, this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc)

	num_features = data.shape[1]

	print("data: ", data.shape)
	print("data.shape: ", data.shape)

	def getPredictions(scores, edges, clf):
		out = []
		if clf.classifier_select == "RF":
			pred_prob = clf.predict_proba(scores)
			pred_class = clf.predict(scores)
		elif clf.classifier_select == "CNN":
			num_samples, num_scores, num_frc = scores.shape
			# print("data: ", data.shape)
			scores = scores.reshape((num_samples, num_scores, num_frc, 1))
			pred_prob = clf.predict_proba(scores)
			pred_class = clf.predict(scores)
		elif clf.classifier_select == "LS":
			num_samples, num_scores, num_frc = scores.shape
			# print("data: ", data.shape)
			scores = scores.reshape((num_samples, num_scores*num_frc))
			pred_class = clf.clf.predict(scores)
			pred_class = pred_class.reshape((len(pred_class), 1))
			pred_prob = clf.clf.predict_proba(scores)
			pred_prob = pred_prob[np.arange(len(pred_prob)), (pred_prob[:,0] < pred_prob[:, 1]).astype(int)]
			pred_prob = pred_prob.reshape((len(pred_prob), 1))
		print("pred_class: ", len(pred_class))
		for i, prediction in enumerate(pred_class):
			# if prediction == 1:
			# print("prediction: ", prediction)
			out.append("%s\t%f" % (edges[i], pred_prob[i]))	#Alternative code that also print label:out.append("%s\t%f\t%i" % (edges[i], pred_prob[i], prediction))
		print("out: ", len(out))
		return out

	out = []
	if clf.classifier_select == "RF":
		tmpscores = np.zeros((100000, num_features))
	elif clf.classifier_select == "CNN" or clf.classifier_select == "LS":
		tmpscores = np.zeros((100000, num_ep, num_frc))
	edges = [""]*100000
	k = 0
	chunk_num=1
	scoreCalc.open()
	print "to predict: %i" % scoreCalc.to_predict
	for line in range(scoreCalc.to_predict):
		if k % 100000==0 and k != 0:
			out.extend(getPredictions(tmpscores[0:k, :], edges[0:k], clf))

			if clf.classifier_select == "RF":
				tmpscores = np.zeros((100000, num_features))
			elif clf.classifier_select == "CNN" or clf.classifier_select == "LS":
				tmpscores = np.zeros((100000, num_ep, num_frc))
			edges = [""] * 100000
			if verbose:
				print "Completed chunk %i" % chunk_num
				chunk_num += 1
			k = 0
		edge, edge_scores = scoreCalc.get_next()
		if edge == "" or edge_scores == []: continue
		# print("edge_scores: ", edge_scores)
		# edge_scores = edge_scores.reshape(1, -1)
		if clf.classifier_select == "RF":
			edge_scores = edge_scores.reshape(1, -1)
		edges[k] = edge
		tmpscores[k,0:(edge_scores.shape)[1]] = edge_scores
		k += 1
	scoreCalc.close()
	out.extend(getPredictions(tmpscores[0:k,:], edges[0:k], clf))
	return out

def get_FA_data(anno_source, taxid, file="", datadir = ""):
	functionalData = ""
	if anno_source == "GM":

		genemania = CS.Genemania(taxid)
		#genemania = CS.Genemania("6239")
		functionalData = genemania.getScoreCalc()

	elif anno_source == "STRING":

		string = CS.STRING(taxid, datadir)
		functionalData = string.getScoreCalc()

	elif anno_source == "FILE":
		if file == "":
			print "When using FILE tag please suppy path to file containing functional annotation using -F file+path"
			sys.exit()
		# the supplied functional evidence data needs to have the correct header row...
		externaldata = CS.ExternalEvidence(file)
		#externaldata.readFile()
		functionalData = externaldata.getScoreCalc()

	else:
		print "EPIC only support GeneMane, STRING, and flat file input please use the followign tags for anno_source GM, STRING, FILE. Returning empty string object."
	return functionalData

def make_predictions(score_calc, mode, clf, gs, output_dir, fun_anno="", verbose = False, folds = 5, train_test_ratio = None, num_ep = 2, num_frc = 27):
	mode = mode.upper()

	def get_edges_from_network(network):
		out = {}
		for edge in network:
			edge, score = edge.rsplit("\t", 1)
			out[edge] = score
		return out

	networks = []
	# predicts using experiment only
	if mode == "EXP" or mode == "BR": networks.append(predictInteractions(score_calc, clf, gs, output_dir, True, verbose, folds, train_test_ratio, num_ep, num_frc))
	#predicts using fun_anno only
	if mode == "FA"or mode == "BR":
		if fun_anno=="":
			# TODO make illigal argument error
			print "if using only functional annotation for prediction functional annotation (fun_anno param != "") must not be empty"
			sys.exit()
		networks.append(predictInteractions(score_calc, clf, gs, output_dir, True, verbose, folds, train_test_ratio, num_ep, num_frc))

	#predict using both functional annotation and exp
	if mode == "COMB" or mode == "BR":
		tmp_score_calc = copy.deepcopy(score_calc)
		print tmp_score_calc.scores.shape
		tmp_score_calc.add_fun_anno(fun_anno)
		print tmp_score_calc.scores.shape
		networks.append(predictInteractions(score_calc, clf, gs, output_dir, True, verbose, folds, train_test_ratio, num_ep, num_frc))

	# return error when no networks is predicted
	if len(networks) == 0:
		print "Error no networks predicted"
		sys.exit()
	# return finised network when only one network is predicted, which happens in any mode expect final
	elif len(networks) ==1:
		return networks[0]
	# use bias reduced method to merge experimental, functional annotation, and combined network
	else:
		exp = get_edges_from_network(networks[0])
		fa = get_edges_from_network(networks[1])
		merged = get_edges_from_network(networks[2])
		br_edges = set(exp.keys()) | (set(merged.keys()) - set(fa.keys()))
		br_network = []
		network_edges = []
		for edge in br_edges:
			if edge in exp: score = exp[edge]
			else: score = merged[edge]
			br_network.append("%s\t%s" % (edge, score))
			network_edges.append(edge)
		return br_network

# a fucntion added by Lucas HU, for testing algorithm stability
# just return the edges of PPIs without interaction scores
# a trial version
def get_network_edges(network):
	network_edges = []

	for items in network:
		proteinA, proteinB, score = items.split("\t")
		edge = "\t".join(sorted([proteinA, proteinB]))
		network_edges.append(edge)

	return network_edges


def predict_clusters(predF, outF):
	dir_path = os.path.dirname(os.path.realpath(__file__))
	clustering_CMD = "java -jar %s/cluster_one-1.0.jar %s > %s" % (dir_path, predF, outF)
	os.system(clustering_CMD)












def load_data(data, scores, orthmap="", fc=2, mfc=1):

	if type(data) is list:
		paths = data
	else:
		paths = [os.path.join(data,fn) for fn in next(os.walk(data))[2]]

	elutionDatas = []
	elutionProts = set([])
	for elutionFile in paths:
		if elutionFile.rsplit(os.sep, 1)[-1].startswith("."): continue
		elutionFile = elutionFile.rstrip()
		elutionData = CS.ElutionData(elutionFile, frac_count=fc, max_frac_count=mfc)
		if orthmap !="":
			if orthmap != False:
				mapper = GS.Inparanoid("", inparanoid_cutoff=1)
				mapper.readTable(orthmap, direction=0)
				elutionData.orthmap(mapper)
		elutionDatas.append(elutionData)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())
		for score in scores:
			score.init(elutionData)
	return elutionProts, elutionDatas


def get_reference_from_net(target_taxid):
	if target_taxid != "9606":
		reference_clusters = [GS.Intact_clusters(True), GS.CORUM(True), GS.QuickGO("9606", True), GS.QuickGO(target_taxid, False)]
	else:
		reference_clusters = [GS.Intact_clusters(False), GS.CORUM(False), GS.QuickGO("9606", False)]
	return reference_clusters

def create_goldstandard(clusters, target_taxid, valprots, pos_neg_ratio):
	if target_taxid !="9606" and target_taxid != "":
		orthmap = GS.Inparanoid(taxid=target_taxid)
	else:
		orthmap = ""

	gs = GS.Goldstandard_from_Complexes("Goldstandard", pos_neg_ratio)
	gs.make_reference_data(clusters, orthmap, found_prots=valprots)
	return gs


def clustering_evaluation(eval_comp, pred_comp, prefix, verbose= True):
	head = "\t".join(["%s%s" % (prefix, h) for h in ["mmr", "overlapp", "simcoe", "mean_simcoe_overlap", "sensetivity", "ppv", "accuracy", "sep"]])

	if len(pred_comp.complexes) > 0:
		cluster_scores = "\t".join(map(str, pred_comp.clus_eval(eval_comp)))
	else:
		cluster_scores =  "\t".join(["0"]*8)
	composite_score = 0
	if verbose:
		tmp_head = head.split("\t")
		tmp_scores = cluster_scores.split("\t")
		#composite_score = 0
		for i in range(len(tmp_head)):
			print "%s\t%s" % (tmp_head[i], tmp_scores[i])

			# add composite score output.
			# added by Lucas HU, a trial function.
			if tmp_head[i] == "mmr" or tmp_head[i] == "overlapp" or tmp_head[i] == "accuracy":
				composite_score = composite_score + float(tmp_scores[i])

		print "composite score is: " + str(composite_score)

	return cluster_scores, head, composite_score

def clusters_to_json(clusters, network, frac_names, eData):
	graph = {}
	for line in network:
		edge, score = line.rsplit("\t", 1)
		graph[edge] = score

	cy_elements = []
	nodes = [] # used to send entwork to cytoscape
	edges = [] # used to send entwork to cytoscape

	net_nodes = set([])
	for complex in clusters.complexes:
		prots = list(clusters.complexes[complex])
		for i in range(len(prots)):
			protA = prots[i]
			for j in range(i+1, len(prots)):
				protB = prots[j]
				edge = "\t".join(sorted([protA, protB]))
				score = 0.5
				if edge in graph: score = graph[edge]
				nodeA = "%s_%s" % (protA, str(complex))
				nodeB = "%s_%s" % (protB, str(complex))
				net_nodes.add(nodeA)
				net_nodes.add(nodeB)
				edge = {
					'group': 'edges',
					'data': {
						'source': nodeA,
						'target': nodeB,
						'score': float(score),
					}
				}
				cy_elements.append(edge)
				edges.append(edge)

	for gene in net_nodes:
		name, cluster_id = gene.split("_")
		node = {
			'group': 'nodes',
			'data': {
				'id': str(gene),
				'name': name,
				'cluster_id': cluster_id
			}}
		for i in range(len(frac_names)):
			score = 0
			if name in eData: score = eData[name][i]
			node['data'][frac_names[i]] = float(score)
		cy_elements.append(node)
		nodes.append(node)

	return json.dumps(cy_elements, default=lambda cy_elements: cy_elements.__dict__), edges, nodes


def json_to_cy_js(div_id, json_str):
	return """

	                $('#cy').show();
	                var cy = window.cy = cytoscape({
	                    container: document.getElementById('%s'),
	                    layout: { },
	                    elements: %s,
	                    style: [
	                     {
	                        selector: 'node',
	                        style: {
	                          'content': 'data(name)',
	                          'font-size': 12,
	                          'text-valign': 'center',
	                          'text-halign': 'center',
	                          'background-color': '#555',
	                          'text-outline-color': '#555',
	                          'text-outline-width': 1.75,
	                          'color': '#fff',
	                          'overlay-padding': 6,
	                          'z-index': 10
	                        }
	                      },
	                      {
	                        selector: 'edge',
	                        style: {
                              'content': 'data(score)',
	                          'line-color': 'black',
	                          'width': 'mapData(score, 0.5, 1, 0, 20)',
	                        }
	                      },
	                    ]
	                });

                    cy.elements().components().forEach( (eles, i, components) => {
                    let n = Math.floor( Math.sqrt( components.length ) );
                    let w = 2000; // width of bb for 1 cmp
                    let h = 2000; // height "

                    eles.makeLayout({
                        name: 'circle',
                        boundingBox: {
                        x1: w * (i %% n),
      x2: w * (1 + (i %% n)),  // this line fixed
      y1: Math.floor(i / n) * h,
      y2: (1 + Math.floor(i / n)) * h
                    }
                    }).run();
                    });
	            """ % (div_id, json_str)

def elutionDatas_to_treeview(eDatas, foundprots, normed=False):
	out = {}
	colnums = {}
	header = []
	all_prots = set([])
	for eData in eDatas:
		name = eData.name
		colnum = eData.elutionMat.shape[1]
		colnums[name] = colnum
		prefix = "%s.F" % name
		header.extend(map(lambda x: "%s%s" % (prefix, x), range(1,colnum+1)))
		all_prots |= set(eData.prot2Index.keys())


	if foundprots != "": all_prots &= foundprots

	for prot in all_prots:
		out[prot] = []
		for eData in eDatas:
			scores = [0]*colnums[eData.name]
			if eData.hasProt(prot):
				scores = eData.getElution(prot, normed)
			out[prot].extend(scores)
	return header, out


def prep_network_for_cy(nodes, edges):
	# Basic Setup
	PORT_NUMBER = 1234
	# IP = '192.168.1.1'
	IP = 'localhost'
	BASE = 'http://' + IP + ':' + str(PORT_NUMBER) + '/v1/'

	# Header for posting data to the server as JSON
	HEADERS = {'Content-Type': 'application/json'}
	network_cy = {
		'data': {'name': "EPIC clusters"},
		"elements": {"nodes": nodes, 'edges': edges}
	}
	return BASE, json.dumps(network_cy), HEADERS

# a fucntion added by Lucas HU to test the stability of prediction using n-fold cross validation
# focus on the PPI level, and see if each time predicetd similar set of PPIs, we use n_fold of data to do this...
def stability_evaluation(n_fold, all_gs, scoreCalc, clf, output_dir, mode, anno_source, taxid, anno_F):

	tmp_train_eval_container = (all_gs.split_into_n_fold2(n_fold, set(scoreCalc.ppiToIndex.keys()))["turpleKey"])

	#create the dictionary to store the predicted PPIs
	PPIs_dict_for_each_fold = {}

	#create the dictionary to store the predicted complexes
	complexes_dict_for_each_fold = {}

	for index in range(n_fold):

		train, eval = tmp_train_eval_container[index]

		print "All comp:%i" % len(all_gs.complexes.complexes)
		print "Train comp:%i" % len(train.complexes.complexes)
		print "Eval comp:%i" % len(eval.complexes.complexes)

		print "Num valid ppis in training pos: %i" % len(train.positive)
		print "Num valid ppis in training neg: %i" % len(train.negative)
		print "Num valid ppis in eval pos: %i" % len(eval.positive)
		print "Num valid ppis in eval neg: %i" % len(eval.negative)

		# Evaluate classifier
		bench_clf(scoreCalc, train, eval, clf, output_dir, verbose=True)

		functionalData = ""
		if mode != "exp":
			functionalData = get_FA_data(anno_source, taxid, anno_F)
			print functionalData.scores.shape

		print "the functional evidence data shape is: "


		# Predict protein interaction based on n_fold cross validation
		# network = make_predictions(scoreCalc, "exp", clf, train, output_dir, fun_anno="", verbose = False, folds, train_test_ratio)
		network = make_predictions(scoreCalc, "exp", clf, train, output_dir, fun_anno="", verbose = False)

		# need to write the network into a file for later-on complexes prediction.
		outFH = open("%s.%s.pred.txt" % (output_dir, mode + anno_source), "w")
		print >> outFH, "\n".join(network)
		outFH.close()

		PPIs_dict_for_each_fold[index] = set(get_network_edges(network))

		#predicted_clusters from the predicted PPI network
		predict_clusters("%s.%s.pred.txt" % (output_dir, mode + anno_source),
							   "%s.%s.clust.txt" % (output_dir, mode + anno_source))

		pred_clusters = GS.Clusters(False)
		pred_clusters.read_file("%s.%s.clust.txt" % (output_dir, mode + anno_source))

		complexes_dict_for_each_fold[index] = pred_clusters

		print "fold " + str(index+1) + "is done"

	#create a matrix for storing overlapped matrix, each element in the matrix is a zero.
	overlapped_ratio_matrix_PPIs = np.zeros((n_fold,n_fold))
	overlapped_ratio_matrix_complexes =  np.zeros((n_fold,n_fold))

	for i in range(0, n_fold):
		for j in range(0, n_fold):

			overlapped_ratio_matrix_PPIs[i,j] = (len(PPIs_dict_for_each_fold[i] & PPIs_dict_for_each_fold[j])) / ((len(PPIs_dict_for_each_fold[i]) + len(PPIs_dict_for_each_fold[j])) / 2)

			# calculate the overlapped complexes numbers from both direction and then get the avergae of them
			overlapped_no1 = complexes_dict_for_each_fold[i].getOverlapp(complexes_dict_for_each_fold[j], cutoff = 0.25)
			overlapped_no2 = complexes_dict_for_each_fold[j].getOverlapp(complexes_dict_for_each_fold[i], cutoff = 0.25)

			averaged_overlapped_complexes_no = (overlapped_no1 + overlapped_no2) / 2

			overlapped_ratio_matrix_complexes[i,j] = averaged_overlapped_complexes_no / ((len(complexes_dict_for_each_fold[i].get_complexes()) + len(complexes_dict_for_each_fold[j].get_complexes())) / 2)

	print overlapped_ratio_matrix_PPIs
	print overlapped_ratio_matrix_complexes

	# create the txt file to save the overlap matrix for stabilit testing.
	filename1 = output_dir + " n_fold_corss_validation_PPIs overlap matrix.txt"
	filename2 = output_dir + " n_fold_corss_validation_complexes overlap matrix.txt"

	np.savetxt(filename1, overlapped_ratio_matrix_PPIs, delimiter = '\t')
	np.savetxt(filename2, overlapped_ratio_matrix_complexes, delimiter='\t')

# a function will automatically generate all possible correlation combinations
# all possible correlation combinations were written in a list.
# A trial function added by Lucas Hu
def generate_all_corr_combination(n = 8):

	i = np.array(np.indices(n * (2,))).reshape(n, -1)
	i[:, np.argsort(i.sum(0)[::-1], kind='mergesort')].T[::-1]

	features_list = list()
	for index in range(1, 256):

		features = i[:, index]

		selected_list = list()

		for j in range(0, 8):
			selected_list.append(features[j])

		feature_combination = ''.join(str(e) for e in selected_list)

		features_list.append(feature_combination)

	return features_list


def Goldstandard_from_cluster_File(gsF, pos_neg_ratio, foundprots=""):
	clusters = GS.Clusters(need_to_be_mapped=False)
	clusters.read_file(gsF)
	if foundprots != "": clusters.remove_proteins(foundprots)
	gs = GS.Goldstandard_from_Complexes("All", pos_neg_ratio)
	gs.complexes = clusters
	gs.make_pos_neg_ppis()
	return gs


def eval_plotting(outDir, this_targets, preds, probs, precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc):
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	print("outDir: ", outDir)
	plotCurves([("", curve_roc)], outDir + "/roc.pdf", "False Positive rate", "True Positive Rate")
	recall_vals, precision_vals, threshold = curve_pr
	plotCurves([("", (precision_vals, recall_vals))], outDir + "/pr.pdf", "Recall", "Precision")
	rownames = ["Precision", "Recall", "F-Measure", "AUC PR", "AUC ROC"]
	threshold = np.append(threshold, 1)
	plotCurves([("Precision", (precision_vals, threshold)), ("Recall", (recall_vals, threshold))], outDir + "/cutoff.pdf", "Cutoff", "Evaluation metric score")

	print("     ** Plot confusion_matrix: ")
	# fig, ax = plt.subplots(figsize=(30, 30))
	skplt.metrics.plot_confusion_matrix(this_targets, preds, normalize=False, text_fontsize=30)
	plt.savefig(os.path.join(outDir, "confusion_matrix.png"), dpi=400)

	plt.close()

	print("     ** Plot normalized confusion_matrix: ")
	skplt.metrics.plot_confusion_matrix(this_targets, preds, normalize=True, text_fontsize=30)
	plt.savefig(os.path.join(outDir, "confusion_matrix_normalized.png"), dpi=400)
	plt.close()

	print("     ** Plot curve_roc: ")
	# print("probs: ", probs)
	probs = np.array(probs).T
	# print("probs: ", probs)
	probs_concat = np.vstack((1-probs, probs)).T
	# print("probs_concat: ", probs_concat)
	skplt.metrics.plot_roc(this_targets, probs_concat, classes_to_plot = [1], plot_micro=False, plot_macro=False)
	plt.savefig(os.path.join(outDir, "curve_roc.png"), dpi=400)
	plt.close()

	print("     ** Plot curve_pr: ")
	skplt.metrics.plot_precision_recall(this_targets, probs_concat, classes_to_plot = [1], plot_micro=False)
	plt.savefig(os.path.join(outDir, "curve_pr.png"), dpi=400)
	plt.close()
