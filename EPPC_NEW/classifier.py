import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Conv2D, MaxPooling2D, Flatten, Activation
from keras.utils import to_categorical
from keras.callbacks import EarlyStopping

import scikitplot as skplt
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn import metrics

from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support, plot_precision_recall_curve, plot_roc_curve

from sklearn.metrics import recall_score,accuracy_score
from sklearn.metrics import precision_score,f1_score

from sklearn.datasets import make_classification
import os
import copy

from sklearn import metrics
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support

from .evaluation import get_metrics, plotCurves, classifier_bench_marking

def CNN_pcc_table_diagonal_model(num_scores, num_frc):
    model = Sequential()
    # add model layers
    model.add(Conv2D(128, kernel_size=3, activation='relu', input_shape=(num_scores,num_frc,1)))
    # model.add(MaxPooling2D(pool_size = (2, 2)))
    # model.add(Conv2D(32, kernel_size=3, activation='relu'))
    # model.add(MaxPooling2D(pool_size = (2, 2)))
    # model.add(Conv2D(32, kernel_size=2, activation='relu'))
    model.add(Flatten())
    model.add(Dense(20, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    return model

def CNN_pcc_table_model(num_scores, num_frc):
    model = Sequential()
    # add model layers
    model.add(Conv2D(64, kernel_size=3, activation='relu', input_shape=(num_scores,num_frc,1)))
    model.add(MaxPooling2D(pool_size = (2, 2)))
    model.add(Conv2D(32, kernel_size=3, activation='relu'))
    model.add(MaxPooling2D(pool_size = (2, 2)))
    # model.add(Conv2D(32, kernel_size=2, activation='relu'))
    model.add(Flatten())
    # model.add(Dense(20, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    return model

def CNN_raw_ef_model(num_scores, num_frc):
    # model = Sequential()
    # # add model layers
    # model.add(Conv2D(64, kernel_size=2, activation='relu', input_shape=(num_scores,num_frc,1)))
    # # model.add(Conv2D(32, kernel_size=2, activation='relu'))
    # model.add(Flatten())
    # # model.add(Dense(20, activation='relu'))
    # model.add(Dense(10, activation='relu'))
    # model.add(Dense(5, activation='relu'))
    # model.add(Dense(1, activation='sigmoid'))

    model = Sequential()
    # add model layers
    model.add(Conv2D(128, kernel_size=2, activation='relu', input_shape=(num_scores,num_frc,1)))
    # model.add(Conv2D(64, kernel_size=1, activation='relu'))
    model.add(Flatten())
    model.add(Dense(20, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    return model

def CNN_scores_model(num_scores, num_frc):
    model = Sequential()
    # add model layers
    model.add(Conv2D(64, kernel_size=2, activation='relu', input_shape=(num_scores,num_frc,1)))
    # model.add(Conv2D(32, kernel_size=2, activation='relu'))
    model.add(Flatten())
    # model.add(Dense(20, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    return model

def CNN_classifier(x, y, epoch_num, output_dir, output_dir_2, model_type, test_size=0.1):
    print("%%% test_size: ", test_size)
    num_samples, num_scores, num_frc = x.shape
    x = x.reshape((num_samples, num_scores, num_frc, 1))

    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_size, random_state=0)
    print("     x_train: ", x_train.shape)
    print("     x_test: ", x_test.shape)
    print("     y_train: ", y_train.shape)
    print("     y_test: ", y_test.shape)

    if model_type == "CNN_pcc_table_model":
        model = CNN_pcc_table_model(num_scores, num_frc)
    elif model_type == "CNN_raw_ef_model":
        model = CNN_raw_ef_model(num_scores, num_frc)
    elif model_type == "CNN_scores_model":
        model = CNN_scores_model(num_scores, num_frc)
    elif model_type == "CNN_pcc_table_diagonal_model":
        print("CNN_pcc_table_diagonal_model:")
        model = CNN_pcc_table_diagonal_model(num_scores, num_frc)

    ######################
    ### Model Training ###
    ######################
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    early_stopping = EarlyStopping(monitor='val_loss', patience=50, verbose=2)
    folds = 5
    skf = StratifiedKFold(folds)
    probs = []
    probs_whole = []
    preds = []
    this_targets = []
    i = 1
    ## Assign training X and training Y
    data = x_train
    targets = y_train
    for train, test in skf.split(data, targets):
        print("Processing data...")
        i += 1
        model.fit(data[train], targets[train], validation_data=(data[test], targets[test]), epochs=epoch_num, callbacks=[early_stopping])
        y_pred_test_probas = model.predict(data[test])
        y_pred_test_probas_whole = np.vstack(((1-y_pred_test_probas)[:,0], y_pred_test_probas[:,0])).T
        y_pred_test = copy.deepcopy(y_pred_test_probas)
        y_pred_test[y_pred_test >= 0.5] = 1
        y_pred_test[y_pred_test < 0.5] = 0
        probs.extend(y_pred_test_probas)
        probs_whole.extend(y_pred_test_probas_whole)
        preds.extend(y_pred_test)
        this_targets.extend(targets[test])

    ##################################
    ### Cross-validation evalution ###
    ##################################
    precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = get_metrics(probs, preds, this_targets)

    if not os.path.exists(output_dir + "/Cross_Validation"):
        os.makedirs(output_dir + "/Cross_Validation")

    format = "pdf"
    plotCurves([("", curve_roc)], output_dir + "/Cross_Validation/roc." + format, "False Positive rate", "True Positive Rate")
    recall_vals, precision_vals, threshold = curve_pr
    plotCurves([("", (precision_vals, recall_vals))], output_dir + "/Cross_Validation/pr." + format, "Recall", "Precision")
    rownames = ["Precision", "Recall", "F-Measure", "AUC PR", "AUC ROC"]
    threshold = np.append(threshold, 1)
    plotCurves([("Precision", (precision_vals, threshold)), ("Recall", (recall_vals, threshold))], output_dir + "/Cross_Validation/cutoff." + format, "Cutoff", "Evaluation metric score")

    print("     ** Plot confusion_matrix: ")
    skplt.metrics.plot_confusion_matrix(this_targets, preds, normalize=False)
    plt.savefig(os.path.join(output_dir, "Cross_Validation", "confusion_matrix.png"), dpi=400)
    plt.close()

    print("     ** Plot normalized confusion_matrix: ")
    skplt.metrics.plot_confusion_matrix(this_targets, preds, normalize=True)
    plt.savefig(os.path.join(output_dir, "Cross_Validation", "confusion_matrix_normalized.png"), dpi=400)
    plt.close()

    print("     ** Plot curve_roc: ")
    skplt.metrics.plot_roc(this_targets, probs_whole, classes_to_plot = [1], plot_micro=False, plot_macro=False)
    plt.savefig(os.path.join(output_dir, "Cross_Validation", "curve_roc.png"), dpi=400)
    plt.close()

    print("     ** Plot curve_pr: ")
    skplt.metrics.plot_precision_recall(this_targets, probs_whole, classes_to_plot = [1], plot_micro=False)
    plt.savefig(os.path.join(output_dir, "Cross_Validation", "curve_pr.png"), dpi=400)
    plt.close()

    #################################
    ### Validation data evalution ###
    #################################
    y_test_pred_probas = model.predict(x_test)
    y_test_pred_probas_whole = np.vstack(((1-y_test_pred_probas)[:,0], y_test_pred_probas[:,0])).T
    y_test_pred = copy.deepcopy(y_test_pred_probas)
    y_test_pred[y_test_pred >= 0.5] = 1
    y_test_pred[y_test_pred < 0.5] = 0
    precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = get_metrics(y_test_pred_probas, y_test_pred, y_test)

    if not os.path.exists(output_dir + "/Validation"):
        os.makedirs(output_dir + "/Validation")

    format = "pdf"
    plotCurves([("", curve_roc)], output_dir + "/Validation/roc." + format, "False Positive rate", "True Positive Rate")
    recall_vals, precision_vals, threshold = curve_pr
    plotCurves([("", (precision_vals, recall_vals))], output_dir + "/Validation/pr." + format, "Recall", "Precision")
    rownames = ["Precision", "Recall", "F-Measure", "AUC PR", "AUC ROC"]
    threshold = np.append(threshold, 1)
    plotCurves([("Precision", (precision_vals, threshold)), ("Recall", (recall_vals, threshold))], output_dir + "/Validation/cutoff." + format, "Cutoff", "Evaluation metric score")

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



    # model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    # early_stopping = EarlyStopping(monitor='val_loss', patience=50, verbose=2)
    # # model.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=epoch_num)
    # model.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=epoch_num, callbacks=[early_stopping])
    #
    # y_pred_test_probas = model.predict(x_test)
    # y_pred_test = copy.deepcopy(y_pred_test_probas)
    # y_pred_test[y_pred_test >= 0.5] = 1
    # y_pred_test[y_pred_test < 0.5] = 0
    # y_pred_test_probas_whole = np.vstack(((1-y_pred_test_probas)[:,0], y_pred_test_probas[:,0])).T
    # # print(y_pred_test)
    # # print(y_pred_test_probas)
    # # print(y_pred_test_probas_whole)
    # classifier_bench_marking(y_test, y_pred_test, y_pred_test_probas, y_pred_test_probas_whole, output_dir)
    #
    #
    #
    # y_pred_probas = model.predict(x)
    # y_pred = copy.deepcopy(y_pred_probas)
    # y_pred[y_pred >= 0.5] = 1
    # y_pred[y_pred < 0.5] = 0
    # y_pred_probas_whole = np.vstack(((1-y_pred_probas)[:,0], y_pred_probas[:,0])).T
    # # print(y_pred_test)
    # # print(y_pred_test_probas)
    # # print(y_pred_test_probas_whole)
    # classifier_bench_marking(y, y_pred, y_pred_probas, y_pred_probas_whole, output_dir_2)

    return model
















def SVM_classifier(x, y, ratio, output_dir, test_size=0.1):
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_size, random_state=0)

    x_train_num_samples, x_train_num_scores, x_train_num_frc = x_train.shape
    x_test_num_samples, x_test_num_scores, x_test_num_frc = x_test.shape

    x_train = np.reshape(x_train, (x_train_num_samples, x_train_num_scores*x_train_num_frc))
    x_test = np.reshape(x_test, (x_test_num_samples, x_test_num_scores*x_test_num_frc))
    print("     x_train: ", x_train.shape)
    print("     x_test: ", x_test.shape)
    print("     y_train: ", y_train.shape)
    print("     y_test: ", y_test.shape)
    svm = SVC(kernel='linear', probability=True, class_weight={1: 1/ratio})

    svm.fit(x_train, y_train)

    y_pred_test = svm.predict(x_test)
    y_pred_test_probas = svm.predict_proba(x_test)[:,1]
    y_pred_test_probas_whole = svm.predict_proba(x_test)

    # print(len(y_test) - sum(y_pred_test == y_test))
    # print(metrics.accuracy_score(y_test, y_pred_test))

    classifier_bench_marking(y_test, y_pred_test, y_pred_test_probas, y_pred_test_probas_whole, output_dir)
    return svm


def RF_classifier(x, y, ratio, output_dir, output_dir_2, test_size=0.1):
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_size, random_state=0)

    x_train_num_samples, x_train_num_scores, x_train_num_frc = x_train.shape
    x_test_num_samples, x_test_num_scores, x_test_num_frc = x_test.shape

    x_train = np.reshape(x_train, (x_train_num_samples, x_train_num_scores*x_train_num_frc))
    x_test = np.reshape(x_test, (x_test_num_samples, x_test_num_scores*x_test_num_frc))
    print("     x_train: ", x_train.shape)
    print("     x_test: ", x_test.shape)
    print("     y_train: ", y_train.shape)
    print("     y_test: ", y_test.shape)
    # , n_estimators=1000
    # , class_weight={1:100}
    clf = RandomForestClassifier(random_state=0, n_estimators = 1000, n_jobs=4)

    clf.fit(x_train, y_train)

    y_pred_test = clf.predict(x_test)
    y_pred_test_probas = clf.predict_proba(x_test)[:,1]
    y_pred_test_probas_whole = clf.predict_proba(x_test)
    # print("y_pred_test_probas: ", y_pred_test_probas)
    # print("y_pred_test_probas_whole: ", y_pred_test_probas_whole)

    # print(len(y_test) - sum(y_pred_test == y_test))
    # print(metrics.accuracy_score(y_test, y_pred_test))
    classifier_bench_marking(y_test, y_pred_test, y_pred_test_probas, y_pred_test_probas_whole, output_dir)

    x_num_samples, x_num_scores, x_num_frc = x.shape
    x = np.reshape(x, (x_num_samples, x_num_scores*x_num_frc))
    y_pred = clf.predict(x)
    y_pred_probas = clf.predict_proba(x_test)[:,1]
    y_pred_probas_whole = clf.predict_proba(x_test)
    classifier_bench_marking(y, y_pred, y_pred_probas, y_pred_probas_whole, output_dir_2)
    return clf
