import pandas as pd
import pickle
import os


def extract_drugbank_id_from_pubchem(pubchem_cid, ofile):
    """
    This function fetches drugbank_id from PubChem REST API. This function can take several days to finish.

    :param
    pubchem_cid: pubchem compound id
    ofile: output file to save the dictionary containing pubchem_cid and drugbank_id

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
    if len(pubchem_cid) > 0:    # if new compound found
        print("New compounds found, fetch their DrugBank id using PubChem API")
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

    return cid_dbid_dict


def populate_protein_nodes(ifile, ofile):
    """
    Populate protein node fields and create a TSV file

    :param
    ifile: input file generated running the SQL
    ofile: output file containing protein node details

    :return:
    Nothing, it creates ofile.
    """

    pnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [], 'node_synonyms': [],
              'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [], 'units': []}

    # read input file as dataframe and write data to output file
    print("read the input compound activity file")
    i_df = pd.read_csv(ifile, sep="\t")

    print("create a dictionary using compound activity file")
    for col in i_df.columns:
        if col == 'uniprot':
            uniprot = i_df['uniprot'].to_numpy()
            pnodes['node_id'] = ['UNIPROTKB ' + str(v) for v in uniprot]
        elif col == 'protein_symbol':
            pnodes['node_label'] = i_df['protein_symbol'].to_numpy()
        elif col == 'compound_name':
            pnodes['node_definition'] = i_df['compound_name'].to_numpy()
        elif col == 'protein_name':
            pnodes['node_synonyms'] = i_df['protein_name'].to_numpy()
        elif col == 'ensemblid':
            pnodes['node_dbxrefs'] = i_df['ensemblid'].to_numpy()
        else:
            continue

    # these cols are not in the dataframe
    pnodes['node_namespace'] = ['IDG'] * i_df.shape[0]
    pnodes['value'] = ['' for i in range(i_df.shape[0])]
    pnodes['lowerbound'] = ['' for i in range(i_df.shape[0])]
    pnodes['upperbound'] = ['' for i in range(i_df.shape[0])]
    pnodes['units'] = ['' for i in range(i_df.shape[0])]

    # write to an output dataframe
    print('write data to output protein dataframe')
    o_df = pd.DataFrame(pnodes)
    o_df = o_df.drop_duplicates()
    o_df.to_csv(ofile, sep="\t", index=False)


def populate_compound_nodes(ifile, ofile, drugbank_file=None):
    """
    Populate compound node fields and create a TSV file

    :param
    ifile: input file generated running the SQL
    ofile: output file containing protein node details

    :return:
    Nothing, it creates ofile.
    """

    cnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [], 'node_synonyms': [],
              'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [], 'units': []}

    print("fetch drugbank id for pubchem ids")
    i_df = pd.read_csv(ifile, sep="\t")
    if drugbank_file is None:
        drugbank_file = ofile[:ofile.rfind("/")] + "/drugbank_id.pkl"

    pubchemid = i_df['pubchem_cid'].to_numpy(dtype=int)
    dbid_dict = extract_drugbank_id_from_pubchem(pubchemid, drugbank_file)

    print("create a dictionary using compound activity file")
    for col in i_df.columns:
        if col == 'pubchem_cid':
            pubchemid = i_df['pubchem_cid'].to_numpy(dtype=int)
            cnodes['node_id'] = ['PUBCHEM_CID ' + str(v) for v in pubchemid]
        elif col == 'compound_name':
            cnodes['node_label'] = i_df['compound_name'].to_numpy()
        elif col == 'smiles':
            cnodes['node_definition'] = i_df['smiles'].to_numpy()
        else:
            continue

    # these cols are not in the dataframe
    cnodes['node_namespace'] = ['IDG'] * i_df.shape[0]
    cnodes['node_synonyms'] = ['' for i in range(i_df.shape[0])]
    cnodes['node_dbxrefs'] = [dbid_dict[v] for v in pubchemid]
    cnodes['value'] = ['' for i in range(i_df.shape[0])]
    cnodes['lowerbound'] = ['' for i in range(i_df.shape[0])]
    cnodes['upperbound'] = ['' for i in range(i_df.shape[0])]
    cnodes['units'] = ['' for i in range(i_df.shape[0])]

    # write to an output dataframe
    print('write data to output compound dataframe')
    o_df = pd.DataFrame(cnodes)
    o_df = o_df.drop_duplicates()
    o_df.to_csv(ofile, sep="\t", index=False)


# start of the code
# input output files
inpfile = "compound_activity.tsv"
compound_outfile = "compound_nodes.tsv"
protein_outfile = "protein_nodes.tsv"
dbid_file = "drugbank_id.pkl"

# compound nodes
populate_compound_nodes(inpfile, compound_outfile, drugbank_file=dbid_file)
# protein nodes
populate_protein_nodes(inpfile, protein_outfile)
