import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.metrics import matthews_corrcoef, accuracy_score
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


def generate_ml_data(ifile=None):
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
    print(df)

    # all features
    all_features = [p.strip(' ').strip('"') for p in df['features'].unique()]
    print("Total number of features: ", len(all_features))
    features_to_remove = [p.strip(' ').strip('"') for p in df['features_to_remove'].unique()]
    print("Total number of features to remove: ", len(features_to_remove))
    all_features = np.asarray(list(set(all_features).difference(features_to_remove)))
    print("final number of features: ", len(all_features))

    # positive and unlabeled genes
    pos_genes_cui = [p.strip(' ').strip('"') for p in df['positive_genes_cui'].unique()]
    all_genes_cui = [p.strip(' ').strip('"') for p in df['unknown_genes_cui'].unique()]
    # print(len(set(pos_genes_cui).intersection(all_genes_cui)))
    unlab_genes_cui = list(set(all_genes_cui).difference(pos_genes_cui))
    print("positive and unlabeled genes and common genes: ", len(pos_genes_cui), len(unlab_genes_cui), len(all_genes_cui),
          len(set(pos_genes_cui).intersection(unlab_genes_cui)))

    # determine features for each genes
    cui_name_dict = {}
    gene_features_dict = {g: set() for g in all_genes_cui}
    feature_list = list(df['features'])
    gene_list = list(df['unknown_genes_cui'])
    gene_name_list = list(df['unknown_genes_label'])
    print("number of rows to iterate: ", len(feature_list), len(gene_list))

    for i, gene in enumerate(gene_list):
        p = feature_list[i]
        f = p.strip(' ').strip('"')
        g = gene.strip(' ').strip('"')
        g_label = gene_name_list[i].strip('"').strip(' "')
        # print(g_label)
        gene_features_dict[g].add(f)
        cui_name_dict[g] = g_label
        # keep track of progress
        if (i + 1) % 500000 == 0:
            print("number of records processed: ", i + 1)
    print("number of records processed and total genes: ", i + 1, len(gene_features_dict))

    # generate binary features for each gene
    all_gene_list = []
    all_gene_features = np.zeros((len(gene_features_dict), len(all_features)), dtype=int)
    all_gene_labels = []
    j = 0
    for g, f in gene_features_dict.items():
        all_gene_list.append(cui_name_dict[g])
        # all_gene_list.append(g)
        idx = np.isin(all_features, list(f)).nonzero()[0]
        all_gene_features[j][idx] = 1
        if g in pos_genes_cui:
            all_gene_labels.append(1)
            # print(g, len(idx), np.sum(all_gene_features[j]))
        else:
            all_gene_labels.append(0)
        j += 1

    # remove genes without any features
    t = np.sum(all_gene_features, axis=1)
    ix = np.where(t > 0)[0]
    all_gene_features = all_gene_features[ix]
    all_gene_labels = np.asarray(all_gene_labels)[ix]
    all_gene_list = np.asarray(all_gene_list)[ix]
    return all_gene_features, all_gene_labels, all_gene_list, all_features, cui_name_dict


def main():
    """
    Call subroutines
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-inpfile", help="provide input file (outout of CQL)", required=True)
    args = parser.parse_args()
    # inpfile1 = "../../cfde-distillery_data/KG2ML/parkinson_ouptut.csv"

    # generate positive and unknown data
    X, y, genes, features, cui_name_dict = generate_ml_data(ifile=args.inpfile)
    # p_df = pd.DataFrame(X, columns=features)
    # p_df['label'] = y
    # p_df.to_csv("parkinson_ml_data.tsv", sep="\t", index=False)
    # exit()

    X, y, genes = shuffle(X, y, genes, random_state=123)
    print(X.shape, len(y), np.sum(y), len(genes), len(features))

    # ML model
    r = Counter(y)[0] / Counter(y)[1]
    print("Train and test the model (neg/pos): ", r)
    params = {'max_depth': 5, 'n_jobs': 16, 'scale_pos_weight': r, 'random_state': 123, 'eval_metric': 'logloss'}
    bst = XGBClassifier(**params)
    # bst = LogisticRegression(random_state=1234, n_jobs=8, class_weight={0: 1, 1: r})
    # bst = RandomForestClassifier(max_depth=3, random_state=1234, class_weight={0: 1, 1: r}, criterion='log_loss')
    # bst = SVC(gamma='auto', probability=True, class_weight={0: 1, 1: r}, random_state=1234)
    preds = cross_val_predict(bst, X, y, cv=5, method='predict_proba')
    print("AUC-ROC: ", roc_auc_score(y, preds[:, 1]))
    print("MCC: ", matthews_corrcoef(y, np.round(preds[:, 1])))
    print("Accuracy: ", accuracy_score(y, np.round(preds[:, 1])))

    # confusion matrix
    tn, fp, fn, tp = confusion_matrix(y, np.round(preds[:, 1])).ravel()
    print("sensitivity (recall): ", tp / (tp + fn))
    print("precision (PPV): ", tp / (tp + fp))

    # run PU method
    # get parameters from user for PULSNAR algorithm.
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
                                    bin_method='scott', bw_method='hist', lowerbw=0.01, upperbw=0.5, optim='local',
                                    calibration=True, calibration_data='PU', calibration_method='sigmoid',
                                    calibration_n_bins=100, smooth_isotonic=False,
                                    classification_metrics=False,
                                    n_iterations=1, kfold=5, kflips=1,
                                    pulsnar_params_file=user_param_file)

    # get results
    res = pls.pulsnar(X, y, tru_label=y, rec_list=genes)
    print("Estimated alpha: {0}".format(res['estimated_alpha']))


if __name__ == "__main__":
    main()
