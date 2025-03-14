from dpTice import *
import pandas as pd
import argparse
from copy import deepcopy
from sklearn.utils import shuffle
from scipy import stats

def fetch_gene_labels(ifile=None):
    """
    fetch gene CUI and labels from the file. create a dictionary with CUI as key and label as value
    """
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

    # get all positives, unlabeled and features from the dataframe
    pos_genes_cui = [p.strip(' ').strip('"') for p in df['positive_genes_cui']]
    all_genes_cui = [p.strip(' ').strip('"') for p in df['unknown_genes_cui']]
    all_features = [p.strip(' ').strip('"') for p in df['features']]
    features_to_remove = set([p.strip(' ').strip('"') for p in df['features_to_remove']])

    # select features for all genes
    gene_features_dict = {}
    for i, u in enumerate(all_genes_cui):
        p = pos_genes_cui[i]
        f = all_features[i]
        if f not in features_to_remove:
            gene_features_dict.setdefault(p, set()).add(f)
            gene_features_dict.setdefault(u, set()).add(f)

    # determine unique genes and features
    all_gene_list = sorted(list(gene_features_dict.keys()))
    all_feature_list = np.sort(list(set().union(*gene_features_dict.values())))
    print("Unique genes and features count: ", len(all_gene_list), len(all_feature_list))

    # generate binary features for each gene
    all_gene_labels = []
    all_gene_names = []
    all_gene_features = np.zeros((len(all_gene_list), len(all_feature_list)), dtype=np.uint8)

    pos_genes_cui = set(pos_genes_cui)
    for j, g in enumerate(all_gene_list):
        all_gene_names.append(cui_name_dict[g])
        f = list(gene_features_dict[g])
        idx = np.isin(all_feature_list, f).nonzero()[0]
        all_gene_features[j][idx] = 1
        if g in pos_genes_cui:
            all_gene_labels.append(1)
        else:
            all_gene_labels.append(0)

    # remove genes without any features
    gene_feature_count = np.sum(all_gene_features, axis=1)
    # print(min(gene_feature_count[p_idx]), max(gene_feature_count[p_idx]), np.mean(gene_feature_count[p_idx]), np.median(gene_feature_count[p_idx]))
    ix = np.where(gene_feature_count >= round(np.median(gene_feature_count)))[0]  # to ensure gene has at least one feature
    return all_gene_features[ix], np.asarray(all_gene_labels)[ix], np.asarray(all_gene_names)[ix]


def main():
    """
    This code calls KM PU method to estimate alpha.
    """
    # local variables to save values
    n_iterations = 40
    alphas = []

    # get command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-inpfile", help="provide input file (output of CQL)", required=True, default="../../CFDE_data/KG2ML/cerebral_all.csv")
    parser.add_argument("-gene_label_file", help="provide input file containing genes and their labels", required=True,
                        default="../../CFDE_data/KG2ML/gene_label_data.csv")
    args = parser.parse_args()

    # fetch labels for all genes
    orig_cui_name_dict = fetch_gene_labels(ifile=args.gene_label_file)

    # generate positive and unknown data
    orig_X, orig_y, orig_genes = generate_ml_data(orig_cui_name_dict, ifile=args.inpfile)
    print("Data shape: ", orig_X.shape, len(orig_y), np.sum(orig_y), len(orig_genes))

    for itr in range(n_iterations):
        print("\nRunning models for iteration: {0}".format(itr + 1))

        # Classification performance before applying PUL
        X, Y, genes = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_genes)
        X, Y = shuffle(X, Y, random_state=itr)

        # call TICE algorithm for alpha estimation
        tice_alpha = tice_wrapper(X, Y)
        print("TiCE alpha: {0}".format(tice_alpha))
        alphas.append(tice_alpha)

    print("\nmean TiCE alpha and 95% CI: ", np.mean(alphas), stats.t.interval(0.95, len(alphas) - 1, loc=np.mean(alphas), scale=stats.sem(alphas)))


if __name__ == "__main__":
    main()
