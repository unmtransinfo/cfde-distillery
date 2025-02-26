import pandas as pd
import numpy as np
from collections import Counter


def select_unique_sab(ifile):
    """
    Read the input file and select unique SAB
    """
    unique_sab = set()
    df = pd.read_csv(ifile, sep="\t", header=0)
    nodes = df['Node'].to_numpy()

    # select SAB
    for i, node in enumerate(nodes):
        v = node.split(":")[0]
        if 'PUBCHEM_CID' in v:
            v = 'PUBCHEM'
        if 'EFO' in v:
            v = 'EFO'
        unique_sab.add(v)
        if (i + 1) % 500000 == 0:
            print("nodes read and unique SAB: ", i + 1, len(unique_sab))

    print("nodes read and unique SAB: ", i + 1, len(unique_sab))
    return unique_sab


def get_data_from_codes_file(ifile, unique_sab_list):
    """
    Read CODEs.csv and return a dictionary
    """
    df_codes = pd.read_csv(ifile, sep=",", header=0, low_memory=False)
    code_ids = df_codes['CodeID:ID'].to_numpy()
    sabs = df_codes['SAB'].to_numpy()
    codes = df_codes['CODE'].to_numpy()
    values = df_codes['value:float'].to_numpy()
    lowerbounds = df_codes['lowerbound:float'].to_numpy()
    upperbounds = df_codes['upperbound:float'].to_numpy()
    code_units = df_codes['unit'].to_numpy()

    # create a dictionary
    code_vals_dict = {}
    for i, code_id in enumerate(code_ids):
        if sabs[i] in unique_sab_list:
            if code_id not in code_vals_dict:
                code_vals_dict[code_id] = [sabs[i], codes[i], values[i], lowerbounds[i], upperbounds[i], code_units[i]]
            else:
                print("duplicate key found: ", code_id)
        else:
            continue

        # keep track of processing
        if (i + 1) % 500000 == 0:
            print("records read: ", i + 1)
    print("records read: ", i + 1)
    return code_vals_dict


def get_data_from_cui_codes_file(ifile, unique_sab_list):
    """
    Read CUI-CODEs.csv and return a dictionary
    """
    df_cui_codes = pd.read_csv(ifile, sep=",", header=0, low_memory=False)
    cuis = df_cui_codes[':START_ID'].to_numpy()
    codes = df_cui_codes[':END_ID'].to_numpy()

    code_cui_dict = {}
    for i, code in enumerate(codes):
        sab = code.split(":")[0]
        if sab not in unique_sab_list:
            continue
        else:
            if code not in code_cui_dict:
                code_cui_dict[code] = [cuis[i]]
            else:
                code_cui_dict[code].append(cuis[i])

        # keep track of processing
        if (i + 1) % 500000 == 0:
            print("records read: ", i + 1)
    print("records read: ", i + 1)

    return code_cui_dict


def get_data_from_cui_sui_file(ifile):
    """
    Read CUI-SUIs.csv and return a dictionary
    """
    cui_synonyms_dict = {}

    df_cui_sui = pd.read_csv(ifile, sep=",", header=0)
    cuis = df_cui_sui[':START_ID'].to_numpy()
    suis = df_cui_sui[':END_ID'].to_numpy()

    for i, cui in enumerate(cuis):
        if cui not in cui_synonyms_dict:
            cui_synonyms_dict[cui] = [suis[i]]
        else:
            cui_synonyms_dict[cui].append(suis[i])

        # keep track of processing
        if (i + 1) % 500000 == 0:
            print("records read: ", i + 1)
    print("records read: ", i + 1)

    return cui_synonyms_dict


def create_nodes_file(code_vals_dict, code_cui_dict, cui_synonyms_dict, node_file):
    """
    create TSV file containing all nodes
    """
    rec_cnt = 0
    unique_cui_list = set()
    with open(node_file, "w") as fo:
        hdr = "CUI\tSAB\tCode\tLabel\tValue\tLower_bound\tUpper_bound\tUnit\tSynonyms\n"
        fo.write(hdr)
        for code_id, vals in code_vals_dict.items():
            sab, code, val, lb, ub, unit = vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]

            # select CUI and node synonyms
            cuis = code_cui_dict[code_id]
            node_synonyms = []
            for cu in cuis:
                cui = cu
                if cu in cui_synonyms_dict:
                    node_synonyms.extend([v for v in cui_synonyms_dict[cu] if str(v) != 'nan'])
                else:
                    continue
            if len(node_synonyms) > 0:
                label = node_synonyms[0]
                synonym = ' | '.join(node_synonyms)
            else:
                label, synonym = '', ''

            line = str(cui) + "\t" + str(sab) + "\t" + str(code) + "\t" + str(label) + "\t" + str(val) + "\t" + str(lb) + "\t" + str(ub) + "\t" + str(
                unit) + "\t" + str(synonym) + "\n"
            unique_cui_list.add(cui)
            line = line.replace('nan', ' ')
            fo.write(line)
            rec_cnt += 1

            # keep track of processing
            if rec_cnt % 250000 == 0:
                print("nodes written: ", rec_cnt)
    print("nodes written: ", rec_cnt)
    return unique_cui_list


def create_edges_file(code_cui_dict, ifile, edge_file, all_cuis):
    """
    Read CUI-CUIs.csv and create edge file
    """
    df_cui_cui = pd.read_csv(ifile, sep=",", header=0, low_memory=False)
    source_cuis = df_cui_cui[':START_ID'].to_numpy()
    target_cuis = df_cui_cui[':END_ID'].to_numpy()
    relations = df_cui_cui[':TYPE'].to_numpy()
    sabs = df_cui_cui['SAB'].to_numpy()
    evidences = df_cui_cui['evidence_class:string'].to_numpy()

    # convert code_cui_dict to cui_code_dict
    cui_code_dict = {}
    for code, cuis in code_cui_dict.items():
        for cui in cuis:
            cui_code_dict[cui] = code

    rec_cnt = 0
    print("write edge data to a file")
    with open(edge_file, "w") as fo:
        hdr = "source\trelation\ttarget\tsource_label\ttarget_label\tevidence\tSAB\n"
        fo.write(hdr)
        for i, su in enumerate(source_cuis):
            if su in all_cuis and target_cuis[i] in all_cuis:
                if su in cui_code_dict:
                    source_label = cui_code_dict[su]
                else:
                    source_label = ''

                if target_cuis[i] in cui_code_dict:
                    target_label = cui_code_dict[target_cuis[i]]
                else:
                    target_label = ''

                line = str(su) + "\t" + str(relations[i]) + "\t" + str(target_cuis[i]) + "\t" + str(source_label) + "\t" + str(
                    target_label) + "\t" + str(evidences[i]) + "\t" + str(sabs[i]) + "\n"
                line = line.replace('nan', ' ')
                fo.write(line)
                rec_cnt += 1

                # keep track of processing
                if rec_cnt % 500000 == 0:
                    print("edges written: ", rec_cnt)
            else:
                continue
    print("edges written: ", rec_cnt)


def separate_node_files(unique_sabs, ifile):
    """
    split nodes files based on SAB
    """
    df_nodes = pd.read_csv(ifile, sep="\t", header=0, low_memory=False)
    for sab in unique_sabs:
        df1 = df_nodes[df_nodes['SAB'] == sab]
        ofile = "../results/" + sab + ".tsv"
        df1.to_csv(ofile, index=False, sep="\t")
        print("{0} file written".format(ofile))


def main():
    """
    main function
    """
    base_input = "../data/neo4j/import/"
    # select unique SAB
    inpfile = "../results/all_nodes.tsv"
    unique_sab_list = select_unique_sab(inpfile)
    unique_sab_list.update(['IDGP', 'IDGD', 'LINCS', 'GO', 'NCBI'])
    print(unique_sab_list)
    exit()
    # read files
    # file 1
    print("\nRead CODEs.csv and process records")
    code_vals_dict = get_data_from_codes_file(base_input + "CODEs.csv", unique_sab_list)
    print("number of records in code_vals_dict: ", len(code_vals_dict))

    # file 2
    print("\nRead CUI-CODEs.csv and process records")
    code_cui_dict = get_data_from_cui_codes_file(base_input + "CUI-CODEs.csv", unique_sab_list)
    print("number of records in code_cui_dict: ", len(code_cui_dict))

    # file 3
    print("\nRead CUI-SUIs.csv and process records")
    cui_synonyms_dict = get_data_from_cui_sui_file(base_input + "CUI-SUIs.csv")
    print("number of records in cui_synonyms_dict: ", len(cui_synonyms_dict))

    # create node file
    print("\ngenerate nodes.tsv")
    node_file = "../results/nodes_pk.tsv"
    all_cuis_list = create_nodes_file(code_vals_dict, code_cui_dict, cui_synonyms_dict, node_file)
    print("number of unique CUIs in all nodes: ", len(all_cuis_list))

    # create edge file
    # file 4
    print("\nRead CUI-CUIs.csv and create edge file")
    inpfile = base_input + "CUI-CUIs.csv"
    edge_file = "../results/edges_pk.tsv"
    create_edges_file(code_cui_dict, inpfile, edge_file, all_cuis_list)

    # split nodes into multiple files
    print("\nseparate nodes by SAB")
    nodes_inpfile = "../results/nodes_pk.tsv"
    separate_node_files(unique_sab_list, nodes_inpfile)


if __name__ == "__main__":
    main()
