import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_predict
from collections import Counter
from sklearn.utils import shuffle
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score, matthews_corrcoef, accuracy_score, brier_score_loss, f1_score, average_precision_score
from scipy.sparse import csr_matrix
from PULSNAR import PULSNAR
import os
import argparse
from copy import deepcopy
import yaml

class KG2MLEstimator:
    """
    the KG2ML pipeline utilizes PU learning method to discover novel associations among biomedical entities
    """
    def __init__(self, param_file=None, is_csr=True, is_scar=True):
        if param_file is not None:  # this is parameter file for PULSNAR package
            self.user_param_file = param_file
        else:
            self.user_param_file = "testparams/pu_alpha_params.yaml"
        self.is_csr = is_csr
        self.is_scar = is_scar

    def fetch_xgboost_parameters(self):
        """
        from the PULSNAR parameter file, fetch all parameters for the XGBoost
        """
        with open(self.user_param_file, 'r') as f:
            pu_params = yaml.safe_load(f)
        return pu_params['XGB_params']

    def fetch_biomed_entity_cui_name(self, ifile=None):
        """
        fetch CUI and names from the file. create a dictionary with CUI as key and name as value
        """
        print(f"fetch CUI and names reading input file: {ifile}")
        # read the input file
        if ifile is not None:
            df = pd.read_csv(ifile)
        else:
            raise ValueError("input file is missing")
        df.columns = df.columns.str.replace(' ', '')

        cui_label_dict = {v.strip().strip('"'): df['gene_label'][i].strip().strip('"') if df['gene_label'][i].strip().strip('"') != 'NULL' else
        v.strip().strip('"').split(' CUI')[0] for i, v in enumerate(df['gene_CUI'])}

        return cui_label_dict

    def generate_ml_data(self, cui_name_dict, ifile=None, min_feature_count=5):
        """
        read input file and generate data for PU method
        """
        print(f"generate ML data using input file: {ifile}")
        # all biomedical entities
        if ifile is not None:
            df = pd.read_csv(ifile)
        else:
            raise ValueError("input file is missing")
        df.columns = df.columns.str.replace(' ', '')

        # get all positives, unlabeled and features from the dataframe
        pos_biomed_entities_cui = [p.strip(' ').strip('"') for p in df['positive_genes_cui']]
        all_biomed_entities_cui = [p.strip(' ').strip('"') for p in df['unknown_genes_cui']]
        all_features = [p.strip(' ').strip('"') for p in df['features']]
        features_to_remove = set([p.strip(' ').strip('"') for p in df['features_to_remove']])

        # select features for all entities
        biomed_entity_features_dict = {}
        for i, u in enumerate(all_biomed_entities_cui):
            p = pos_biomed_entities_cui[i]
            f = all_features[i]
            if f not in features_to_remove:
                biomed_entity_features_dict.setdefault(p, set()).add(f)
                biomed_entity_features_dict.setdefault(u, set()).add(f)

        # determine unique biomed entity and features
        all_biomed_entity_list = sorted(list(biomed_entity_features_dict.keys()))
        all_feature_list = np.sort(list(set().union(*biomed_entity_features_dict.values())))
        print("Unique biomed entity and features count: ", len(all_biomed_entity_list), len(all_feature_list))

        # generate binary features for each biomed entity
        all_biomed_entity_labels = []
        all_biomed_entity_names = []
        all_biomed_entity_features = np.zeros((len(all_biomed_entity_list), len(all_feature_list)), dtype=np.uint8)

        pos_biomed_entities_cui = set(pos_biomed_entities_cui)
        for j, g in enumerate(all_biomed_entity_list):
            all_biomed_entity_names.append(cui_name_dict[g])
            f = list(biomed_entity_features_dict[g])
            idx = np.isin(all_feature_list, f).nonzero()[0]
            all_biomed_entity_features[j][idx] = 1
            if g in pos_biomed_entities_cui:
                all_biomed_entity_labels.append(1)
            else:
                all_biomed_entity_labels.append(0)

        # remove biomed entities without minimum features
        biomed_entity_feature_count = np.sum(all_biomed_entity_features, axis=1)
        # ix = np.where(biomed_entity_feature_count >= round(np.median(biomed_entity_feature_count)))[0]
        ix = np.where(biomed_entity_feature_count >= min_feature_count)[0]  # to ensure biomed entity has at least minimum number of features
        return csr_matrix(all_biomed_entity_features[ix]), np.asarray(all_biomed_entity_labels)[ix], np.asarray(all_biomed_entity_names)[ix], all_feature_list


    def run_pu_classifier(self, X, y, recs, rseed=0, fl_obj=None, v_iter=1):
        """
        Run PU method to estimate alpha
        """
        print(f"\nrun PULSCAR algorithm to estimate alpha value")
        # check if results folder exist. if not, create it
        if not os.path.exists("results"):
            os.makedirs("results")

        # keep original copy of the data
        orig_X, orig_y, orig_recs = deepcopy(X), deepcopy(y), deepcopy(recs)

        # instantiate PULSNAR Classifier
        pls = PULSNAR.PULSNARClassifier(scar=self.is_csr, csrdata=self.is_csr, classifier='xgboost',
                                        bin_method='rice', bw_method='hist', lowerbw=0.001, upperbw=0.5, optim='local',
                                        calibration=True, calibration_data='U', calibration_method='sigmoid',
                                        calibration_n_bins=100, smooth_isotonic=False,
                                        classification_metrics=False,
                                        n_iterations=1, kfold=5, kflips=1,
                                        pulsnar_params_file=self.user_param_file)

        # get results
        X, y, recs = shuffle(X, y, recs, random_state=rseed)
        res = pls.pulsnar(X, y, tru_label=y, rec_list=recs)
        print(f"for iteration number {v_iter}, estimated alpha: {res['estimated_alpha']}")

        # compute classification performance using original data -- BEFORE RUNNING PULSCAR
        print("BEFORE PU - run XGBoost on the original labels to determine classification performance")
        X, y, recs = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_recs)
        X, y = shuffle(X, y, random_state=rseed)
        r = Counter(y)[0] / Counter(y)[1]
        print(f"BEFORE PU - unlabeled to positive ratio: {r}")
        xgb_params = self.fetch_xgboost_parameters()
        xgb_params['scale_pos_weight'] = r
        xgb_params['random_state'] = rseed
        # print(f"BEFORE PU - parameters for XGBoost: {xgb_params}")
        bst = XGBClassifier(**xgb_params)
        preds = cross_val_predict(bst, X, y, cv=5, method='predict_proba')

        # performance metrics - BEFORE RUNNING PULSCAR
        tn, fp, fn, tp = confusion_matrix(y, np.round(preds[:, 1])).ravel()
        b_auc_val = roc_auc_score(y, preds[:, 1])
        b_mcc_val = matthews_corrcoef(y, np.round(preds[:, 1]))
        b_acc_val = accuracy_score(y, np.round(preds[:, 1]))
        b_brier_val = brier_score_loss(y, preds[:, 1])
        b_f1_val = f1_score(y, np.round(preds[:, 1]))
        b_sensitivity_val = tp / (tp + fn)
        b_specificity_val = tn / (tn + fp)
        b_precision_val = tp / (tp + fp)
        b_aps_val = average_precision_score(y, preds[:, 1])

        print("AFTER PU - run XGBoost on the imputed labels to determine classification performance")
        X, y, recs = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_recs)
        df = pd.read_csv(res['prediction_file'], sep="\t", header=0)
        unlab_rec_ids = df['rec_id'].to_numpy()
        calibrated_prob = df['calibrated_prob'].to_numpy()

        # select top alpha*U unlabeled examples as probable positives
        idx = np.argsort(calibrated_prob)[::-1][:int(res['estimated_alpha'] * len(calibrated_prob))]  # select top alpha*U unlabeled
        probable_positives = unlab_rec_ids[idx]
        print(f"number of probable positives among {len(unlab_rec_ids)} unlabeled records: {len(probable_positives)}")

        # flip labels for probable positives
        idx = np.isin(recs, probable_positives).nonzero()[0]
        y[idx] = 1
        print(f"number of true + probable positives: {len(np.where(y==1)[0])}, number of probable negatives: {len(np.where(y==0)[0])}")

        # ML model on imputed labels
        X, y = shuffle(X, y, random_state=rseed)
        r = Counter(y)[0] / Counter(y)[1]
        print(f"AFTER PU - unlabeled to positive ratio: {r}")
        xgb_params['scale_pos_weight'] = r
        # print(f"AFTER PU - parameters for XGBoost: {xgb_params}")
        bst = XGBClassifier(**xgb_params)
        preds = cross_val_predict(bst, X, y, cv=5, method='predict_proba')

        # performance metrics - BEFORE RUNNING PULSCAR
        tn, fp, fn, tp = confusion_matrix(y, np.round(preds[:, 1])).ravel()
        a_auc_val = roc_auc_score(y, preds[:, 1])
        a_mcc_val = matthews_corrcoef(y, np.round(preds[:, 1]))
        a_acc_val = accuracy_score(y, np.round(preds[:, 1]))
        a_brier_val = brier_score_loss(y, preds[:, 1])
        a_f1_val = f1_score(y, np.round(preds[:, 1]))
        a_sensitivity_val = tp / (tp + fn)
        a_specificity_val = tn / (tn + fp)
        a_precision_val = tp / (tp + fp)
        a_aps_val = average_precision_score(y, preds[:, 1])

        # write performance metrics to a file
        line_out = (str(v_iter) + "\t" + str(b_auc_val) + "\t" + str(a_auc_val) + "\t" + str(b_mcc_val) + "\t" + str(a_mcc_val) + "\t" +
                    str(b_acc_val) + "\t" + str(a_acc_val) + "\t" + str(b_brier_val) + "\t" + str(a_brier_val) + "\t" +
                    str(b_f1_val) + "\t" + str(a_f1_val) + "\t" + str(b_sensitivity_val) + "\t" + str(a_sensitivity_val) + "\t" +
                    str(b_specificity_val) + "\t" + str(a_specificity_val) + "\t" + str(b_precision_val) + "\t" + str(a_precision_val) + "\t" +
                    str(b_aps_val) + "\t" + str(a_aps_val) + "\t" + str(res['estimated_alpha']) + "\n")
        fl_obj.write(line_out)

        # get important features
        model = XGBClassifier(**xgb_params).fit(X, y)
        feature_importance = model.get_booster().get_score(importance_type='gain')

        # print classification performance metrics
        print(f"AUC - before PU: {b_auc_val}, after PU: {a_auc_val}")
        print(f"MCC - before PU: {b_mcc_val}, after PU: {a_mcc_val}")
        print(f"Accuracy - before PU: {b_acc_val}, after PU: {a_acc_val}")
        print(f"Brier Loss - before PU: {b_brier_val}, after PU: {a_brier_val}")
        print(f"Sensitivity - before PU: {b_sensitivity_val}, after PU: {a_sensitivity_val}")
        print(f"Specificity - before PU: {b_specificity_val}, after PU: {a_specificity_val}")
        print(f"Precision - before PU: {b_precision_val}, after PU: {a_precision_val}")
        print(f"APS - before PU: {b_aps_val}, after PU: {a_aps_val}")
        return res, feature_importance

def main():
    """
    This code calls PU classifier to estimate the proportion of positives among unlabeled examples; determines those probable positives
    using their calibrated probabilities. It also computes performance metrics before and after identifying probable positives.
    """
    min_feature_count = 1  # discard a record if it does have enough features

    # get command line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-n_iterations", type=int, default=40)
    parser.add_argument("-iofiles", help="provide yaml file containing io files", default="io_data.yaml")
    parser.add_argument("-pu_params_file", help="provide yaml file containing PULSNAR parameters", default="testparams/pu_alpha_params.yaml")
    p_args = parser.parse_args()
    with open(p_args.iofiles, 'r') as fi:
        iodata = yaml.safe_load(fi)

    # instantiate KG2ML estimator
    kg2mle = KG2MLEstimator(param_file=p_args.pu_params_file, is_csr=True, is_scar=True)

    # fetch names for all selected biomedical entities
    orig_cui_name_dict = kg2mle.fetch_biomed_entity_cui_name(ifile=iodata['names_data'])
    print(f"number of elements in orig_cui_name_dict: {len(orig_cui_name_dict)}")

    # generate positive and unknown data
    orig_X, orig_y, orig_recs, orig_features = kg2mle.generate_ml_data(orig_cui_name_dict, ifile=iodata['input_data'], min_feature_count=min_feature_count)
    print(f"Data shape: {orig_X.shape}, {len(orig_y)}, positive count: {np.where(orig_y == 1)[0].shape[0]}, unlabeled count: {np.where(orig_y == 0)[0].shape[0]}")

    # start processing
    dir_path = os.path.dirname(iodata['cl_performance_file'])   # get directory path of the file
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)

    alphas,rec_preds, imp_features_gain = [], {}, {}
    with open(iodata['cl_performance_file'], 'w') as fout:
        file_hdr = "iteration\tbefore_AUC\tafter_AUC\tbefore_MCC\tafter_MCC\tbefore_Accuracy\tafter_Accuracy\tbefore_Brier_Loss\tafter_Brier_Loss\tbefore_F1\tafter_F1\tbefore_Sensitivity\tafter_Sensitivity\tbefore_Specificity\tafter_Specificity\tbefore_Precision\tafter_Precision\tbefore_APS\tafter_APS\tAlpha\n"
        fout.write(file_hdr)
        for itr in range(p_args.n_iterations):
            X, y, recs = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_recs)
            pu_res, feature_importance_gain = kg2mle.run_pu_classifier(X, y, recs, rseed=1234*(itr+1), fl_obj=fout, v_iter=itr+1)

            # fetch alpha and predictions
            alphas.append(pu_res['estimated_alpha'])
            df = pd.read_csv(pu_res['prediction_file'], sep="\t", header=0)
            pred_recs = df['rec_id'].to_numpy()
            calibrated_prob = df['calibrated_prob'].to_numpy()
            predicted_prob = df['predicted_prob'].to_numpy()

            for j, rec in enumerate(pred_recs):
                rec_preds.setdefault(rec, {})
                rec_preds[rec].setdefault('calibrated_prob', []).append(calibrated_prob[j])
                rec_preds[rec].setdefault('predicted_prob', []).append(predicted_prob[j])

            # fetch important features and gain scores
            for kk, vv in feature_importance_gain.items():
                imp_features_gain.setdefault(kk, []).append(vv)

    # save mean predictions
    print(f"mean estimated alpha: {np.mean(alphas)}")
    with open(iodata['pu_prediction_file'], 'w') as fout:
        line_hdr = "rec_name\tpredicted_prob\tcalibrated_prob\n"
        fout.write(line_hdr)
        for b_rec, vals in rec_preds.items():
            line_out = str(b_rec) + "\t" + str(np.mean(vals['predicted_prob'])) + "\t" + str(np.mean(vals['calibrated_prob'])) + "\n"
            fout.write(line_out)

    # save mean gain scores for each of important features
    with open(iodata['pu_feature_importance_file'], 'w') as fout:
        line_hdr = "feature_name\tgain_value\n"
        fout.write(line_hdr)
        for f_rec, vals in imp_features_gain.items():
            j = int(f_rec[1:])  # XGBoost returns f1, f2, ... as features if names are not provided
            line_out = str(orig_features[j]) + "\t" + str(np.mean(vals)) + "\n"
            fout.write(line_out)


if __name__ == "__main__":
    main()
