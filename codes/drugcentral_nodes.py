import pandas as pd


def populate_drugcentral_nodes(ifile, ofile):
    """
    Populate node fields and create a TSV file
    """

    pnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [], 'node_synonyms': [],
              'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [], 'units': []}

    # read input file as dataframe and write data to output file
    print("read the input compound activity file")
    i_df = pd.read_csv(ifile, sep="\t")

    print("create a dictionary using compound activity file")
    for col in i_df.columns:
        if col == 'snomed_conceptid':
            snomed_id = i_df['snomed_conceptid'].to_numpy()
            pnodes['node_id'] = ['SNOMEDCT_US ' + str(v) for v in snomed_id]
        elif col == 'omop_concept_name':
            pnodes['node_label'] = i_df['omop_concept_name'].to_numpy()
        elif col == 'omop_concept_id':
            omop_concept = i_df['omop_concept_id'].to_numpy()
            pnodes['node_dbxrefs'] = ['OMOP ' + str(v) for v in omop_concept]
        else:
            continue

    # these cols are not in the dataframe
    pnodes['node_namespace'] = ['IDG'] * i_df.shape[0]
    pnodes['node_definition'] = ['' for i in range(i_df.shape[0])]
    pnodes['node_synonyms'] = ['' for i in range(i_df.shape[0])]
    pnodes['value'] = ['' for i in range(i_df.shape[0])]
    pnodes['lowerbound'] = ['' for i in range(i_df.shape[0])]
    pnodes['upperbound'] = ['' for i in range(i_df.shape[0])]
    pnodes['units'] = ['' for i in range(i_df.shape[0])]

    # write to an output dataframe
    print('write data to output protein dataframe')
    o_df = pd.DataFrame(pnodes)
    o_df = o_df.drop_duplicates()
    o_df.to_csv(ofile, sep="\t", index=False)


# start of the code
inpfile = "/home/praveen/Documents/work/cfde-distillery_data/drug_disease_input.tsv"
outfile = "/home/praveen/Documents/work/cfde-distillery_data/nodes_disease.tsv"
populate_drugcentral_nodes(inpfile, outfile)
