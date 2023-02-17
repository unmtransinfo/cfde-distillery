import pandas as pd


def populate_edges_data(ifile, ofile):
    """
    Populate edge fields and create a TSV file
    """

    # dict to store values
    edges = {'subject_id': [], 'relationship': [], 'object_id': [], 'evidence_class': []}

    # read input file as dataframe and write data to output file
    print("read the input compound activity file")
    i_df = pd.read_csv(ifile, sep="\t")

    print("create a dictionary using compound activity file")
    for col in i_df.columns:
        if col == 'umls_cui':
            uniprot = i_df['umls_cui'].to_numpy()
            edges['object_id'] = ['UMLS ' + str(v) for v in uniprot]
        elif col == 'pubchem_cid':
            pubchemid = i_df['pubchem_cid'].to_numpy(dtype=int)
            edges['subject_id'] = ['PUBCHEM_CID ' + str(v) for v in pubchemid]
        else:
            continue

    edges['relationship'] = ['indication'] * i_df.shape[0]
    edges['evidence_class'] = ['' for i in range(i_df.shape[0])]

    # write to an output dataframe
    print('write data to output edge dataframe')
    o_df = pd.DataFrame(edges)
    o_df = o_df.drop_duplicates()
    o_df.to_csv(ofile, sep="\t", index=False)


# start of the code
inpfile = "/home/praveen/Documents/work/cfde-distillery_data/drug_disease.csv"
outfile = "/home/praveen/Documents/work/cfde-distillery_data/indication_edges.tsv"
populate_edges_data(inpfile, outfile)