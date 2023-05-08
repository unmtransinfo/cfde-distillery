from DistilleryNodesEdges import DistilleryNodesEdges
import numpy as np
import pandas as pd


def validate_nodes_edges_data(p_nodes_file, c_nodes_file, d_nodes_file, cd_edge_file, cp_edge_file):
    """
    Validate nodes and edges
    """
    MISSING_DATA = False
    # protein nodes
    p_df = pd.read_csv(p_nodes_file, sep="\t", header=0)
    n_uniprot_ids = set(p_df['node_id'].to_numpy())

    # compound nodes
    c_df = pd.read_csv(c_nodes_file, sep="\t", header=0)
    n_pubchem_ids = set(c_df['node_id'].to_numpy())

    # disease nodes
    d_df = pd.read_csv(d_nodes_file, sep="\t", header=0)
    n_snomed_ids = set(d_df['node_id'].to_numpy())

    # compound disease edges
    cd_df = pd.read_csv(cd_edge_file, sep="\t", header=0)
    e_snomed_ids = set(cd_df['object_id'].to_numpy())
    e1_pubchem_ids = set(cd_df['subject_id'].to_numpy())

    # compound protein edges
    cp_df = pd.read_csv(cp_edge_file, sep="\t", header=0)
    e_uniprot_ids = set(cp_df['object_id'].to_numpy())
    e2_pubchem_ids = set(cp_df['subject_id'].to_numpy())

    # validate compound disease edges
    if len(e_snomed_ids.difference(n_snomed_ids)) > 0:
        print("Validation failed!!! compound disease edge file has additional SNOMED ids")
        MISSING_DATA = True
    if len(e1_pubchem_ids.difference(n_pubchem_ids)) > 0:
        print("Validation failed!!! compound disease edge file has additional PubChem ids")
        MISSING_DATA = True

    # validate compound protein edges
    if len(e_uniprot_ids.difference(n_uniprot_ids)) > 0:
        print("Validation failed!!! compound protein edge file has additional Uniprot ids")
        MISSING_DATA = True
    if len(e2_pubchem_ids.difference(n_pubchem_ids)) > 0:
        print("Validation failed!!! compound protein edge file has additional PubChem ids")
        MISSING_DATA = True

    return MISSING_DATA


def main():
    """
    Call methods in DistilleryNodesEdges class to generate TSV files for nodes and edges
    """
    # inputs
    tcrd_inpfile = "/home/praveen/Documents/work/cfde-distillery_data/compound_activity_input.tsv"
    dc_inpfile = "/home/praveen/Documents/work/cfde-distillery_data/drug_disease_input.tsv"
    dbid_inpfile = "/home/praveen/Documents/work/cfde-distillery_data/compound_name_drugbank_id.pkl"

    # outputs
    p_nodes_file = "/home/praveen/Documents/work/cfde-distillery_data/nodes_protein.tsv"
    c_nodes_file = "/home/praveen/Documents/work/cfde-distillery_data/nodes_compound.tsv"
    d_nodes_file = "/home/praveen/Documents/work/cfde-distillery_data/nodes_disease.tsv"
    cd_edge_file = "/home/praveen/Documents/work/cfde-distillery_data/edges_indication_compound_disease.tsv"
    cp_edge_file = "/home/praveen/Documents/work/cfde-distillery_data/edges_bioactivity_compound_protein.tsv"

    # instantiate class
    dne = DistilleryNodesEdges(tcrd_ifile=tcrd_inpfile, dc_ifile=dc_inpfile, dbid_names_ifile=dbid_inpfile,
                               compound_node_file=c_nodes_file, protein_node_ofile=p_nodes_file,
                               disease_node_ofile=d_nodes_file, compound_disease_edge_file=cd_edge_file,
                               compound_protein_edge_file=cp_edge_file)

    # generate nodes
    dne.generate_nodes(node_type="protein")
    dne.generate_nodes(node_type="compound")
    dne.generate_nodes(node_type="disease")

    # generate edges
    dne.generate_edges(edge_type="compound_disease")
    dne.generate_edges(edge_type="compound_protein")

    # *********************************************************************************#
    # RUN THE FOLLOWING CODE FOR VALIDATION:                                           #
    # If a node is present (uniprot, pubchem, snomed) in the edge file, it MUST be     #
    # present in the node file(s)                                                      #
    # *********************************************************************************#
    print("Validate data in edges and nodes files")
    missing_info = validate_nodes_edges_data(p_nodes_file, c_nodes_file, d_nodes_file, cd_edge_file, cp_edge_file)
    if missing_info:
        print("FAIL - Nodes files do not have all nodes")
    else:
        print("Validation passed")


if __name__ == "__main__":
    main()
