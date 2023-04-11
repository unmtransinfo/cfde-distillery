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
        if col == 'uniprot':
            uniprot = i_df['uniprot'].to_numpy()
            edges['object_id'] = ['UNIPROTKB ' + str(v) for v in uniprot]
        elif col == 'pubchem_cid':
            pubchemid = i_df['pubchem_cid'].to_numpy(dtype=int)
            edges['subject_id'] = ['PUBCHEM_CID ' + str(v) for v in pubchemid]
        elif col == 'act_type':
            edges['evidence_class'] = i_df['act_type'].to_numpy()
        else:
            continue

    edges['relationship'] = ['bioactivity'] * i_df.shape[0]

    # write to an output dataframe
    print('write data to output edge dataframe')
    o_df = pd.DataFrame(edges)
    o_df = o_df.drop_duplicates()
    o_df.to_csv(ofile, sep="\t", index=False)


# start of the code
inpfile = "/home/praveen/Documents/work/cfde-distillery_data/compound_activity_input.tsv"
outfile = "/home/praveen/Documents/work/cfde-distillery_data/edges_bioactivity_compound_protein.tsv"
populate_edges_data(inpfile, outfile)
