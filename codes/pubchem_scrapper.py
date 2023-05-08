import pandas as pd
import pickle
import os


def fetch_all_pubchem_ids(tfile, dcfile):
    """
    Read both TCRD and DrugCentral files to get all pubchem ids
    """
    i_df = pd.read_csv(tfile, sep="\t")
    pubchemids = set(i_df['pubchem_cid'].to_numpy(dtype=int))
    print("fetched pubchem ids from TCRD compound activity file: ", len(pubchemids))

    i_df = pd.read_csv(dcfile, sep="\t")
    pubchemids.update(set(i_df['pubchem_cid'].to_numpy(dtype=int)))
    print("fetched pubchem ids from DrugCentral disease file: ", len(pubchemids))
    return pubchemids


def extract_drugbank_id_and_name(pubchem_cid, ofile):
    """
    This function fetches drugbank id and all synonyms for a given compound from pubchem REST API.
    This function can take several hours to finish.
    """

    cid_dbid_name_dict = {}
    rcount = 0

    # check if drugbank file exists
    if os.path.exists(ofile):
        with open(ofile, "rb") as fh:
            cid_dbid_name_dict = pickle.load(fh)
        print("Number of elements in pubchem_cid set (original): ", len(pubchem_cid))
        pubchem_cid = set(pubchem_cid).difference(cid_dbid_name_dict.keys())
        print("Number of elements in pubchem_cid set (after subtracting): ", len(pubchem_cid))

    # fetch data from pubchem
    for pbcid in set(pubchem_cid):
        if pbcid not in cid_dbid_name_dict:
            cid_dbid_name_dict[pbcid] = {}

        pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(pbcid) + "/JSON/"
        print(pubchem_url)
        FOUND_DRUG_BANK_ID = False
        FOUND_COMPOUND_NAMES = False

        # read the JSON file
        df = pd.read_json(pubchem_url, orient='index', encoding='latin-1')

        # find name
        for rec in df['Section'][0]:
            if FOUND_COMPOUND_NAMES:
                break
            for k1, v1 in rec.items():
                if FOUND_COMPOUND_NAMES:
                    break
                if k1 == 'TOCHeading' and v1 == 'Names and Identifiers':
                    for vals in rec['Section']:
                        if FOUND_COMPOUND_NAMES:
                            break
                        for k2, v2 in vals.items():
                            if FOUND_COMPOUND_NAMES:
                                break
                            if k2 == 'TOCHeading' and v2 == 'Synonyms':
                                for val in vals['Section']:
                                    if FOUND_COMPOUND_NAMES:
                                        break
                                    for k3, v3 in val.items():
                                        if FOUND_COMPOUND_NAMES:
                                            break
                                        if k3 == 'TOCHeading' and v3 == 'Depositor-Supplied Synonyms':
                                            c_names = "|".join(synm['String'] for synm in val['Information'][0]['Value']['StringWithMarkup'])
                                            cid_dbid_name_dict[pbcid]['names'] = c_names
                                            FOUND_COMPOUND_NAMES = True
        # if name not found, write ''
        if not FOUND_COMPOUND_NAMES:
            cid_dbid_name_dict[pbcid]['names'] = ''

        # find drugbank id
        for rec in df['Reference'][0]:
            if FOUND_DRUG_BANK_ID:
                break
            for k, v in rec.items():
                if FOUND_DRUG_BANK_ID:
                    break
                if k == 'SourceName' and v == 'DrugBank':
                    cid_dbid_name_dict[pbcid]['DrugBank'] = rec['SourceID']
                    FOUND_DRUG_BANK_ID = True

        # if not drug_bank_id, write ''
        if not FOUND_DRUG_BANK_ID:
            cid_dbid_name_dict[pbcid]['DrugBank'] = ''

        rcount += 1
        # data dictionary in a pickle file after every 250 records.
        if rcount % 250 == 0:
            print("saved {0} records in the pickle file".format(rcount))
            with open(ofile, "wb") as fh:
                pickle.dump(cid_dbid_name_dict, fh, protocol=5)

    # save the remaining records
    print("saved {0} records in the pickle file".format(rcount))
    with open(ofile, "wb") as fh:
        pickle.dump(cid_dbid_name_dict, fh, protocol=5)


if __name__ == "__main__":
    # input files
    tcrd_file = "/home/praveen/Documents/work/cfde-distillery_data/compound_activity_input.tsv"
    drugcentral_file = "/home/praveen/Documents/work/cfde-distillery_data/drug_disease_input.tsv"
    # output files
    outfile = "/home/praveen/Documents/work/cfde-distillery_data/compound_name_drugbank_id.pkl"

    print("fetch drugbank id and compound names from pubchem")
    pubchem_ids = fetch_all_pubchem_ids(tcrd_file, drugcentral_file)
    extract_drugbank_id_and_name(pubchem_ids, outfile)
