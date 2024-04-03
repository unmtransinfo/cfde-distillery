import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.metrics import matthews_corrcoef, accuracy_score, brier_score_loss, f1_score, average_precision_score
from collections import Counter
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from scipy.sparse import csr_matrix
from PULSNAR import PULSNAR
import sys
import os
import argparse
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy import stats


def fetch_gene_labels(ifile=None):
    """
    fetch gene CUI and labels from the file. create a dictionary with CUI as key and label as value
    """
    cui_label_dict = {}
    # read the input file
    if ifile is not None:
        df = pd.read_csv(ifile)
    else:
        print("input file1 is missing")
        exit(-1)
    df.columns = df.columns.str.replace(' ', '')

    cui_label_dict = {v.strip().strip('"'): df['gene_label'][i].strip().strip('"') if df['gene_label'][i].strip().strip('"') != 'NULL' else
    v.strip().strip('"').split(' CUI')[0] for i, v in enumerate(df['gene_CUI'])}

    return cui_label_dict


def generate_ml_data(cui_name_dict, ifile=None):
    """
    read input file and generate ML data
    """
    # all genes
    if ifile is not None:
        df = pd.read_csv(ifile)
    else:
        print("input file1 is missing")
        exit(-1)
    df.columns = df.columns.str.replace(' ', '')
    # print(df)

    # all features
    all_features = list(p.strip(' ').strip('"') for p in df['features'].unique())
    print("Total number of features: ", len(all_features))
    features_to_remove = list(p.strip(' ').strip('"') for p in df['features_to_remove'].unique())
    print("Total number of features to remove: ", len(features_to_remove))
    all_features = np.asarray(list(set(all_features).difference(features_to_remove)))
    print("final number of features: ", len(all_features))

    # positive and unlabeled genes
    pos_genes_cui = [p.strip(' ').strip('"') for p in df['positive_genes_cui'].unique()]
    all_genes_cui = [p.strip(' ').strip('"') for p in df['unknown_genes_cui'].unique()]
    # print(len(set(pos_genes_cui).intersection(all_genes_cui)))
    unlab_genes_cui = list(set(all_genes_cui).difference(pos_genes_cui))
    print("positive, unlabeled, and common genes: ", len(pos_genes_cui), len(unlab_genes_cui), len(set(pos_genes_cui).intersection(unlab_genes_cui)))

    # determine features for each genes
    gene_features_dict = {g: set() for g in pos_genes_cui}
    print("gene_features_dict size 1: ", len(gene_features_dict))
    gene_features_dict.update({g: set() for g in unlab_genes_cui})
    print("gene_features_dict size 2: ", len(gene_features_dict))

    feature_list = list(df['features'])
    p_gene_list = list(df['positive_genes_cui'])
    u_gene_list = list(df['unknown_genes_cui'])
    print("number of rows to iterate: ", len(feature_list), len(p_gene_list), len(u_gene_list))

    for i, u_gene in enumerate(u_gene_list):
        p = p_gene_list[i].strip(' ').strip('"')
        u = u_gene.strip(' ').strip('"')
        f = feature_list[i].strip(' ').strip('"')
        gene_features_dict[p].add(f)
        gene_features_dict[u].add(f)

        # keep track of progress
        if (i + 1) % 500000 == 0:
            print("number of records processed: ", i + 1)
    print("number of records processed and total genes: ", i + 1, len(gene_features_dict))

    # generate binary features for each gene
    all_gene_list = []
    all_gene_labels = []
    all_gene_features = np.zeros((len(gene_features_dict), len(all_features)), dtype=int)
    gene_features_dict = dict(sorted(gene_features_dict.items()))  # to ensure the data sequence remains unchanged

    j = 0
    for g, f in gene_features_dict.items():
        all_gene_list.append(cui_name_dict[g])
        idx = np.isin(all_features, list(f)).nonzero()[0]
        all_gene_features[j][idx] = 1
        if g in pos_genes_cui:
            all_gene_labels.append(1)
        else:
            all_gene_labels.append(0)
        j += 1

    # remove genes without any features
    t = np.sum(all_gene_features, axis=1)
    ix = np.where(t > 0)[0]  # to ensure gene has at least one feature
    return all_gene_features[ix], np.asarray(all_gene_labels)[ix], np.asarray(all_gene_list)[ix]


def ml_performance_before_pul(X, y, performance_metrics, params, balanced=False, rseed=0):
    """
    Run XGBoost and compute performance
    """
    disp_msg = "before PU"
    if balanced:
        pos_idx = np.where(y == 1)[0]
        neg_idx = np.where(y == 0)[0]
        np.random.seed(rseed)
        np.random.shuffle(neg_idx)
        sel_neg_idx = neg_idx[:len(pos_idx)]
        idx = np.concatenate([pos_idx, sel_neg_idx])
        print("number of records in the balanced set: ", len(idx))
        # select records for the balanced set
        X = X[idx]
        y = y[idx]
        disp_msg = "[balanced]"

    # ML model
    X, y = shuffle(X, y, random_state=rseed)
    r = Counter(y)[0] / Counter(y)[1]
    # print("Train and test the model (neg/pos): ", r)
    params['scale_pos_weight'] = r
    params['random_state'] = rseed
    bst = XGBClassifier(**params)
    preds = cross_val_predict(bst, X, y, cv=5, method='predict_proba')

    # performance metrics
    auc_val = roc_auc_score(y, preds[:, 1])
    mcc_val = matthews_corrcoef(y, np.round(preds[:, 1]))
    acc_val = accuracy_score(y, np.round(preds[:, 1]))
    brier_val = brier_score_loss(y, preds[:, 1])
    f1_val = f1_score(y, np.round(preds[:, 1]))
    aps_val = average_precision_score(y, preds[:, 1])
    tn, fp, fn, tp = confusion_matrix(y, np.round(preds[:, 1])).ravel()
    recall_val = tp / (tp + fn)
    ppv_val = tp / (tp + fp)

    # update dictionary
    performance_metrics['AUC'].append(auc_val)
    performance_metrics['MCC'].append(mcc_val)
    performance_metrics['Accuracy'].append(acc_val)
    performance_metrics['Brier_loss'].append(brier_val)
    performance_metrics['F1_score'].append(f1_val)
    performance_metrics['APS'].append(aps_val)
    performance_metrics['Sensitivity'].append(recall_val)
    performance_metrics['Precision'].append(ppv_val)

    # display results
    print(f"\nAUC-ROC {disp_msg}: ", auc_val)
    print(f"MCC {disp_msg}: ", mcc_val)
    print(f"Accuracy {disp_msg}: ", acc_val)
    print(f"Brier score loss {disp_msg}: ", brier_val)
    print(f"F1 score {disp_msg}: ", f1_val)
    print(f"APS {disp_msg}: ", aps_val)
    print(f"sensitivity (recall) {disp_msg}: ", recall_val)
    print(f"precision (PPV) {disp_msg}: ", ppv_val)

    return performance_metrics


def run_pu_classifier(X, y, genes, alpha_list, rseed=0):
    """
    Run PU method to estimate alpha
    """
    # update pulsnar_args.yaml if you want to override the default parameters
    if len(sys.argv) < 2:
        user_param_file = 'testparams/gene_alpha.yaml'
    else:
        user_param_file = sys.argv[1]

    # check if results folder exist. if not, create it
    if not os.path.exists("results"):
        os.makedirs("results")

    # instantiate PULSNARClassifier
    pls = PULSNAR.PULSNARClassifier(scar=True, csrdata=False, classifier='xgboost',
                                    bin_method='scott', bw_method='hist', lowerbw=0.1, upperbw=0.5, optim='local',
                                    calibration=True, calibration_data='U', calibration_method='sigmoid',
                                    calibration_n_bins=100, smooth_isotonic=False,
                                    classification_metrics=False,
                                    n_iterations=1, kfold=5, kflips=1,
                                    pulsnar_params_file=user_param_file)

    # get results
    X, y, genes = shuffle(X, y, genes, random_state=rseed)
    res = pls.pulsnar(X, y, tru_label=y, rec_list=genes)
    print("\nEstimated alpha: {0}".format(res['estimated_alpha']))
    alpha_list.append(res['estimated_alpha'])
    return alpha_list


def ml_performance_after_pul(X, y, genes, performance_metrics, params, rseed=0):
    """
    Run XGBoost and compute performance
    """
    disp_msg = "after PU"
    # find probable positives among unlabeled examples
    df = pd.read_csv("predictions.tsv", sep="\t", header=0)
    pred_genes = df['rec_id'].to_numpy()
    pred_calibrated_prob = df['calibrated_prob'].to_numpy()
    idx = np.where(pred_calibrated_prob > 0.5)[0]
    probable_pos_genes = pred_genes[idx]
    print("number of probable positives: ", len(probable_pos_genes))

    # flip labels for probable positives
    idx = np.isin(genes, probable_pos_genes).nonzero()[0]
    print("number of true positives: ", np.sum(y))
    y[idx] = 1
    print("number of true + probable positives: ", np.sum(y))

    # ML model
    X, y, genes = shuffle(X, y, genes, random_state=rseed)
    r = Counter(y)[0] / Counter(y)[1]
    # print("Train and test the model (neg/pos): ", r)
    params['scale_pos_weight'] = r
    params['random_state'] = rseed
    bst = XGBClassifier(**params)
    preds = cross_val_predict(bst, X, y, cv=5, method='predict_proba')

    # performance metrics
    auc_val = roc_auc_score(y, preds[:, 1])
    mcc_val = matthews_corrcoef(y, np.round(preds[:, 1]))
    acc_val = accuracy_score(y, np.round(preds[:, 1]))
    brier_val = brier_score_loss(y, preds[:, 1])
    f1_val = f1_score(y, np.round(preds[:, 1]))
    aps_val = average_precision_score(y, preds[:, 1])
    tn, fp, fn, tp = confusion_matrix(y, np.round(preds[:, 1])).ravel()
    recall_val = tp / (tp + fn)
    ppv_val = tp / (tp + fp)

    # update dictionary
    performance_metrics['AUC'].append(auc_val)
    performance_metrics['MCC'].append(mcc_val)
    performance_metrics['Accuracy'].append(acc_val)
    performance_metrics['Brier_loss'].append(brier_val)
    performance_metrics['F1_score'].append(f1_val)
    performance_metrics['APS'].append(aps_val)
    performance_metrics['Sensitivity'].append(recall_val)
    performance_metrics['Precision'].append(ppv_val)

    # display results
    print(f"\nAUC-ROC {disp_msg}: ", auc_val)
    print(f"MCC {disp_msg}: ", mcc_val)
    print(f"Accuracy {disp_msg}: ", acc_val)
    print(f"Brier score loss {disp_msg}: ", brier_val)
    print(f"F1 score {disp_msg}: ", f1_val)
    print(f"APS {disp_msg}: ", aps_val)
    print(f"sensitivity (recall) {disp_msg}: ", recall_val)
    print(f"precision (PPV) {disp_msg}: ", ppv_val)

    return performance_metrics


def generate_bar_chart(balanced_metrics, before_metrics, after_metrics, metric_list, plt_file="metrics_comparison.png"):
    """
    generate barcharts using the performance metrics
    """
    # Plotting
    fig, ax = plt.subplots()

    # Set positions for the bars
    x = np.arange(len(metric_list))
    bar_width = 0.25

    # Plot bars with error bars
    # rects1 = ax.bar(x, calculate_mean_performance(balanced_metrics, metric_list), bar_width, label='Balanced',
    #                yerr=calculate_95ci(balanced_metrics, metric_list), capsize=2)
    rects1 = ax.bar(x, calculate_mean_performance(before_metrics, metric_list), bar_width, label='XGBoost only',
                    yerr=calculate_95ci(before_metrics, metric_list), capsize=2)
    rects2 = ax.bar(x + bar_width, calculate_mean_performance(after_metrics, metric_list), bar_width, label='XGBoost+PULSCAR',
                    yerr=calculate_95ci(after_metrics, metric_list), capsize=2)

    # Add labels, title, and legend
    ax.set_xlabel('ML model performance metrics', fontsize=12)
    ax.set_ylabel('Scores', fontsize=12)
    ax.set_title('Comparison of metrics with/without PULSCAR', fontsize=13, fontweight='bold')
    ax.set_xticks(x + bar_width / 2)
    ax.set_xticklabels(metric_list, rotation=90)
    ax.legend()
    # ax.grid(visible=True, which='major', axis='y', color='0.25', linestyle='-')
    # ax.grid(visible=True, which='minor', axis='both', color='0.45', linestyle='--')
    ax.minorticks_on()
    ax.xaxis.set_tick_params(which='minor', bottom=False)

    # Save plot as a high-quality PNG file
    plt.savefig(plt_file, dpi=300, bbox_inches='tight')

    # Show plot
    plt.show()


def calculate_mean_performance(performance_metrics, metric_list):
    """
    calculate mean for each of the performance metric
    """
    mean_values = []
    for m in metric_list:
        mean_values.append(np.mean(performance_metrics[m]))
    return mean_values


def calculate_95ci(performance_metrics, metric_list):
    """
    calculate 95% CI for each of the performance metric
    """
    ci_values = []
    for m in metric_list:
        data = performance_metrics[m]
        ci_95p = stats.t.interval(0.95, len(data) - 1, loc=np.mean(data), scale=stats.sem(data))
        ci_values.append([np.mean(data) - ci_95p[0], ci_95p[1] - np.mean(data)])
    return np.transpose(ci_values)


def main():
    """
    This code calls PU classifier to estimate the alpha. Then using true positives, probable positives, and probable negatives, it runs XGBoost
    model and calculates classification performance metrics. It also calculates classification performance metrics without running the PU classifier.
    """
    # loca variables to save values
    n_iterations = 40
    alphas = []
    balanced_metrics = {'AUC': [], 'MCC': [], 'Accuracy': [], 'Brier_loss': [], 'F1_score': [], 'APS': [], 'Sensitivity': [], 'Precision': []}
    before_metrics = {'AUC': [], 'MCC': [], 'Accuracy': [], 'Brier_loss': [], 'F1_score': [], 'APS': [], 'Sensitivity': [], 'Precision': []}
    after_metrics = {'AUC': [], 'MCC': [], 'Accuracy': [], 'Brier_loss': [], 'F1_score': [], 'APS': [], 'Sensitivity': [], 'Precision': []}
    metric_list = ['AUC', 'MCC', 'Accuracy', 'Brier_loss', 'F1_score', 'APS', 'Sensitivity', 'Precision']
    xgb_params = {'max_depth': 5, 'n_jobs': 16, 'eval_metric': 'logloss'}

    # get command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-inpfile", help="provide input file (outout of CQL)", required=True)
    parser.add_argument("-gene_label_file", help="provide input file containing genes and their labels", required=True)
    parser.add_argument("-pltfile", help="provide filename to save bar plot", required=True)
    args = parser.parse_args()

    # fetch labels for all genes
    orig_cui_name_dict = fetch_gene_labels(ifile=args.gene_label_file)

    # generate positive and unknown data
    orig_X, orig_y, orig_genes = generate_ml_data(orig_cui_name_dict, ifile=args.inpfile)
    print("Data shape: ", orig_X.shape, len(orig_y), np.sum(orig_y), len(orig_genes))

    for itr in range(n_iterations):
        print("\nRunning models for iteration: {0}".format(itr + 1))

        # Classification performance before applying PUL
        X, y, genes = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_genes)
        # balanced_metrics = ml_performance_before_pul(X, y, balanced_metrics, xgb_params, balanced=True, rseed=101+itr)  # using balanced set
        before_metrics = ml_performance_before_pul(X, y, before_metrics, xgb_params, balanced=False, rseed=101 + itr)  # all data

        # run PU method
        alphas = run_pu_classifier(X, y, genes, alphas, rseed=101 + itr)

        # classification performance after applying PUL
        X, y, genes = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_genes)
        after_metrics = ml_performance_after_pul(X, y, genes, after_metrics, xgb_params, rseed=101 + itr)

    print("\nbalanced_metrics: ", balanced_metrics)
    print("\nbefore metrics: ", before_metrics)
    print("\nafter metrics: ", after_metrics)
    print("\nalphas: ", alphas)

    # generate barchart using performance metrics
    print("\nmean alpha and 95% CI: ", np.mean(alphas), stats.t.interval(0.95, len(alphas) - 1, loc=np.mean(alphas), scale=stats.sem(alphas)))
    generate_bar_chart(balanced_metrics, before_metrics, after_metrics, metric_list, plt_file=args.pltfile)


if __name__ == "__main__":
    main()
