import pandas as pd

# dict to store values
edges = {'subject_id': [], 'relationship': [], 'object_id': [], 'evidence_class': []}

# read input file as dataframe and write data to output file
print("read the input compound activity file")
i_df = pd.read_csv("../data/compound_activity.tsv", sep="\t")

print("create a dictionary using compound activity file")
for col in i_df.columns:
    if col == 'uniprot':
        uniprot = i_df['uniprot'].to_numpy()
        edges['object_id'] = ['UNIPROTKB ' + str(v) for v in uniprot]
    elif col == 'pubchem_cid':
        edges['subject_id'] = i_df['pubchem_cid'].to_numpy()
    elif col == 'act_type':
        edges['evidence_class'] = i_df['act_type'].to_numpy()
    else:
        continue

edges['relationship'] = ['bioactivity'] * i_df.shape[0]

# write to an output dataframe
print('write data to output edge dataframe')
o_df = pd.DataFrame(edges)
o_df = o_df.drop_duplicates()
o_df.to_csv("../data/edges.tsv", sep="\t", index=False)