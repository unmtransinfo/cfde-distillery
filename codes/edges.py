import pandas as pd
import numpy as np

hdr = "subject_id\trelationship\tobject_id\tevidence_class\n"
nodes_file = open("../data/edges.tsv", "w")
nodes_file.write(hdr)

df = pd.read_csv("../data/compound_activity.tsv", sep="\t")
uniprot_id = df['uniprot'].to_numpy()
act_type = df['act_type'].to_numpy()
pubchem_cid = df['pubchem_cid'].to_numpy()
activity_source = df['activity_source'].to_numpy()

for i in range(df.shape[0]):
    line = str(pubchem_cid[i]) + "\t" + "bioactivity" + "\t" + str("UNIPROTKB ") + str(uniprot_id[i]) + "\t" + \
           str(str(act_type[i])) + "\n"
    nodes_file.write(line)