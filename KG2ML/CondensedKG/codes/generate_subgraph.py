import numpy as np
import pandas as pd
import os


def select_nodes(all_nodes, fnodes):
    """
    select nodes from the file
    """

    for fnode in fnodes:
        if 'PUBCHEM_CID' in fnode:
            v = 'PUBCHEM_CID:' + fnode.strip().split(' ')[1]
            all_nodes.add(v)
        elif 'SNOMEDCT_US' in fnode:
            v = 'SNOMEDCT_US:' + fnode.strip().split(' ')[1]
            all_nodes.add(v)
        elif 'UNIPROTKB' in fnode:
            v = 'UNIPROTKB:' + fnode.strip().split(' ')[1]
            all_nodes.add(v)
        elif 'HGNC' in fnode:
            v = fnode.strip().split(' ')[1]
            all_nodes.add(v)
        elif 'EFO' in fnode:
            all_nodes.add(fnode)
        elif 'GTEXEXP' in fnode:
            v = 'GTEXEXP:' + fnode.strip().split(' ')[1]
            all_nodes.add(v)
        elif 'EXPBINS' in fnode:
            v = 'EXPBINS:' + fnode.strip().split(' ')[1]
            all_nodes.add(v)
        elif 'UBERON' in fnode:
            v = 'UBERON:' + fnode.strip().split(' ')[1]
            all_nodes.add(v)
        else:
            continue

    return all_nodes


# select all nodes
filenames = next(os.walk('../data'))[2]
set_of_nodes = set()
for fl in filenames:
    flname = "../data/" + fl
    dframe = pd.read_csv(flname, sep="\t")

    # select nodes
    print("\nprocessing file: ", flname)
    print("processing columns: ", set(dframe.columns))
    if 'node_id' in dframe.columns:
        nodes = dframe['node_id'].to_numpy()
        set_of_nodes = select_nodes(set_of_nodes, nodes)
    if 'subject_id' in dframe.columns:
        nodes = dframe['subject_id'].to_numpy()
        set_of_nodes = select_nodes(set_of_nodes, nodes)
    if 'object_id' in dframe.columns:
        nodes = dframe['object_id'].to_numpy()
        set_of_nodes = select_nodes(set_of_nodes, nodes)
    if 'uniprot' in dframe.columns:
        nodes = dframe['uniprot'].to_numpy()
        set_of_nodes = select_nodes(set_of_nodes, nodes)
    if 'hgnc' in dframe.columns:
        nodes = dframe['hgnc'].to_numpy()
        set_of_nodes = select_nodes(set_of_nodes, nodes)
    if 'subject' in dframe.columns:
        nodes = dframe['subject'].to_numpy()
        set_of_nodes = select_nodes(set_of_nodes, nodes)
    if 'object' in dframe.columns:
        nodes = dframe['object'].to_numpy()
        set_of_nodes = select_nodes(set_of_nodes, nodes)

# save all nodes
print("save all nodes to a TSV file: ", len(set_of_nodes))
dnodes = pd.DataFrame({'Node': list(set_of_nodes)})
dnodes.to_csv("../results/all_nodes.tsv", sep="\t", index=False)