from sklearn import metrics
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support, plot_precision_recall_curve, plot_roc_curve
import matplotlib.pyplot as plt

from sklearn.metrics import recall_score,accuracy_score
from sklearn.metrics import precision_score,f1_score

from sklearn.datasets import make_classification
import os
import copy

from sklearn import metrics
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support



def get_metrics(probs, preds, targets):
	precision = metrics.precision_score(targets, preds, average=None)[1]
	recall = metrics.recall_score(targets, preds, average=None)[1]
	fmeasure = metrics.f1_score(targets, preds, average=None)[1]
	auc_pr = average_precision_score(targets, preds)
	auc_roc = roc_auc_score(targets, preds)
	curve_pr = precision_recall_curve(targets, probs)
	curve_roc = roc_curve(targets, probs)
	return [precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc]

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


def classifier_bench_marking(y_test, y_pred_test, y_pred_test_probas, y_pred_test_probas_whole, output_dir):
    precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = get_metrics(y_pred_test_probas, y_pred_test, y_test)
    print("     precision: ", precision)
    print("     recall: ", recall)
    print("     fmeasure: ", fmeasure)
    print("     auc_pr: ", auc_pr)
    print("     auc_roc: ", auc_roc)
    print("     curve_pr: ", curve_pr)
    print("     curve_roc: ", curve_roc)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("     ** Plot confusion_matrix: ")
    skplt.metrics.plot_confusion_matrix(y_test, y_test_pred, normalize=False)
    plt.savefig(os.path.join(output_dir, "Validation", "confusion_matrix.png"), dpi=400)
    plt.close()

    print("     ** Plot normalized confusion_matrix: ")
    skplt.metrics.plot_confusion_matrix(y_test, y_test_pred, normalize=True)
    plt.savefig(os.path.join(output_dir, "Validation", "confusion_matrix_normalized.png"), dpi=400)
    plt.close()

    print("     ** Plot curve_roc: ")
    skplt.metrics.plot_roc(y_test, y_test_pred_probas_whole, classes_to_plot = [1], plot_micro=False, plot_macro=False)
    plt.savefig(os.path.join(output_dir, "Validation", "curve_roc.png"), dpi=400)
    plt.close()

    print("     ** Plot curve_pr: ")
    skplt.metrics.plot_precision_recall(y_test, y_test_pred_probas_whole, classes_to_plot = [1], plot_micro=False)
    plt.savefig(os.path.join(output_dir, "Validation", "curve_pr.png"), dpi=400)
    plt.close()
