import pandas as pd
import pickle
import os


def extract_drugbank_id_from_pubchem(pubchem_cid, ofile):
    """
    This function fetches drugbank id from pubchem REST API. This function can take several hours to finish.

    :param
    pubchem_cid: pubchem compound id

    :return:
    a dictionary with pubchem_cid as key and drugbank_id as value
    """
    cid_dbid_dict = {}
    rcount = 0

    # check if drugbank file exists
    if os.path.exists(ofile):
        with open(ofile, "rb") as fh:
            cid_dbid_dict = pickle.load(fh)
        print("Number of elements in pubchem_cid (original): ", len(set(pubchem_cid)))
        pubchem_cid = set(pubchem_cid).difference(cid_dbid_dict.keys())
        print("Number of elements in pubchem_cid (after subtracting): ", len(pubchem_cid))

    # fetch drugbank id
    for pbcid in set(pubchem_cid):
        pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(pbcid) + "/JSON/"
        print(pubchem_url)
        FOUND_DRUG_BANK_ID = False
        # read the JSON file
        df = pd.read_json(pubchem_url, orient='index', encoding='latin-1')
        for rec in df['Reference'][0]:
            for k, v in rec.items():
                if k == 'SourceName' and v == 'DrugBank':
                    cid_dbid_dict[pbcid] = rec['SourceID']
                    FOUND_DRUG_BANK_ID = True
                    break

        # if not drug_bank_id, write ''
        if not FOUND_DRUG_BANK_ID:
            cid_dbid_dict[pbcid] = ''

        rcount += 1
        # data dictionary in a pickle file after every 100 records.
        if rcount % 100 == 0:
            print("saved {0} records in the pickle file".format(rcount))
            with open(ofile, "wb") as fh:
                pickle.dump(cid_dbid_dict, fh, protocol=5)

    # save the remaining records
    print("saved {0} records in the pickle file".format(rcount))
    with open(ofile, "wb") as fh:
        pickle.dump(cid_dbid_dict, fh, protocol=5)


# start of the code
# input output files
inpfile = "/home/praveen/Documents/work/cfde-distillery_data/compound_activity_input.tsv"
dbid_file = "/home/praveen/Documents/work/cfde-distillery_data/drugbank_id.pkl"

print("fetch drugbank id for pubchem ids")
i_df = pd.read_csv(inpfile, sep="\t")
pubchemid = i_df['pubchem_cid'].to_numpy(dtype=int)
extract_drugbank_id_from_pubchem(pubchemid, dbid_file)
