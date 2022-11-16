import pandas as pd

# ******* Protein nodes ******#
pnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [], 'node_synonyms': [],
          'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [], 'units': []}

# read input file as dataframe and write data to output file
print("read the input compound activity file")
i_df = pd.read_csv("../data/compound_activity.tsv", sep="\t")

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
o_df.to_csv("../data/pnodes.tsv", sep="\t", index=False)

# ******* Compound nodes ******#
cnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [], 'node_synonyms': [],
          'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [], 'units': []}

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
cnodes['node_dbxrefs'] = ['' for i in range(i_df.shape[0])]
cnodes['value'] = ['' for i in range(i_df.shape[0])]
cnodes['lowerbound'] = ['' for i in range(i_df.shape[0])]
cnodes['upperbound'] = ['' for i in range(i_df.shape[0])]
cnodes['units'] = ['' for i in range(i_df.shape[0])]

# write to an output dataframe
print('write data to output compound dataframe')
o_df = pd.DataFrame(cnodes)
o_df = o_df.drop_duplicates()
o_df.to_csv("../data/cnodes.tsv", sep="\t", index=False)
