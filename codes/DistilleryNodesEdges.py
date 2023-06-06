import numpy as np
import pandas as pd
import pickle


class DistilleryNodesEdges:
    def __init__(self, tcrd_ifile=None, dc_ifile=None, dbid_names_ifile=None, compound_node_file=None,
                 protein_node_ofile=None, disease_node_ofile=None, compound_disease_edge_file=None,
                 compound_protein_edge_file=None):

        # get arguments
        self.tcrd_ifile = tcrd_ifile  # TCRD file is for compounds and proteins
        self.dc_ifile = dc_ifile  # DrugCentral file is for diseases and compounds
        self.dbid_names_ifile = dbid_names_ifile  # this pickle file should have DrugBank ID and compound names
        self.compound_node_file = compound_node_file  # output file for compound nodes
        self.protein_node_file = protein_node_ofile  # output file for protein nodes
        self.disease_node_file = disease_node_ofile  # output file for disease nodes
        self.compound_disease_edge_file = compound_disease_edge_file  # output file for compound disease edge
        self.compound_protein_edge_file = compound_protein_edge_file  # output file for compound protein edge

    def generate_nodes(self, node_type=None):
        """
        Generate nodes file for the given node type
        """
        if node_type == "protein":
            self.populate_protein_nodes()
        elif node_type == "compound":
            self.populate_compound_nodes()
        elif node_type == "disease":
            self.populate_disease_nodes()
        else:
            print("node type (protein, compound, or disease) is missing")
            exit(-1)

    def generate_edges(self, edge_type=None):
        """
        Generate edges file for the given edge type
        """
        if edge_type == "compound_disease":
            self.populate_compound_disease_edges()
        elif edge_type == "compound_protein":
            self.populate_compound_protein_edges()
        else:
            print("edge type (compound_disease, or compound_protein) is missing")
            exit(-1)

    def populate_protein_nodes(self):
        """
        Populate protein node fields and create a TSV file
        """
        print("\ngenerate protein nodes file")
        pnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [],
                  'node_synonyms': [], 'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [],
                  'units': []}

        if self.tcrd_ifile is None or self.protein_node_file is None:
            print("input/output file is missing")
            exit(-1)

        # read input file as dataframe and write data to output file
        print("fetch data from compound activity (TCRD) file")
        i_df = pd.read_csv(self.tcrd_ifile, sep="\t", header=0)
        i_df = i_df.replace(np.nan, '')  # replace NaN with empty string

        for col in i_df.columns:
            if col == 'uniprot':
                uniprot = list(i_df['uniprot'])
                pnodes['node_id'] = ['UNIPROTKB ' + str(v) for v in uniprot]
            elif col == 'protein_symbol':
                pnodes['node_label'] = list(i_df['protein_symbol'])
            elif col == 'compound_name':
                pnodes['node_definition'] = list(i_df['compound_name'])
            elif col == 'protein_name':
                pnodes['node_synonyms'] = list(i_df['protein_name'])
            elif col == 'ensemblid':
                pnodes['node_dbxrefs'] = list(i_df['ensemblid'])
            else:
                continue

        # these cols are not in the dataframe
        pnodes['node_namespace'] = ['IDG'] * i_df.shape[0]
        pnodes['value'] = ['' for i in range(i_df.shape[0])]
        pnodes['lowerbound'] = ['' for i in range(i_df.shape[0])]
        pnodes['upperbound'] = ['' for i in range(i_df.shape[0])]
        pnodes['units'] = ['' for i in range(i_df.shape[0])]

        # write to an output dataframe
        print('write data to an output TSV file')
        o_df = pd.DataFrame(pnodes)
        o_df = o_df.drop_duplicates()
        o_df.to_csv(self.protein_node_file, sep="\t", index=False)

    def populate_compound_nodes(self):
        """
        Populate compound node fields and create a TSV file
        """
        print("\ngenerate compound nodes file")
        cnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [],
                  'node_synonyms': [], 'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [], 'units': []}

        if self.tcrd_ifile is None or self.dc_ifile is None \
                or self.compound_node_file is None or self.dbid_names_ifile is None:
            print("input/output file is missing")
            exit(-1)

        # fetch dictionary from the pickle file for drugbank ids and names
        with open(self.dbid_names_ifile, 'rb') as fh:
            dbid_dict = pickle.load(fh)

        # read input file as dataframe and write data to output file
        print("fetch data from compound activity (TCRD) file")
        i1_df = pd.read_csv(self.tcrd_ifile, sep="\t", header=0)
        i1_df = i1_df.replace(np.nan, '')  # replace NaN with empty string

        for col in i1_df.columns:
            if col == 'pubchem_cid':
                pubchemid = i1_df['pubchem_cid'].to_numpy(dtype=int)
                cnodes['node_id'] = ['PUBCHEM_CID ' + str(v) for v in pubchemid]
            elif col == 'compound_name':
                pubchemid = i1_df['pubchem_cid'].to_numpy(dtype=int)
                n_labels = []
                for p in pubchemid:
                    if "|" in dbid_dict[p]['names']:
                        n_labels.append(dbid_dict[p]['names'].split("|")[0])
                    else:
                        n_labels.append(dbid_dict[p]['names'])
                cnodes['node_label'] = n_labels
            elif col == 'smiles':
                cnodes['node_definition'] = list(i1_df['smiles'])
            else:
                continue

        # these cols are not in the dataframe
        cnodes['node_namespace'] = ['IDG'] * i1_df.shape[0]
        cnodes['node_synonyms'] = [dbid_dict[p]['names'] for p in pubchemid]
        cnodes['node_dbxrefs'] = [dbid_dict[p]['DrugBank'] for p in pubchemid]
        cnodes['value'] = ['' for i in range(i1_df.shape[0])]
        cnodes['lowerbound'] = ['' for i in range(i1_df.shape[0])]
        cnodes['upperbound'] = ['' for i in range(i1_df.shape[0])]
        cnodes['units'] = ['' for i in range(i1_df.shape[0])]

        # get compounds from DrugCentral file too
        print("fetch data from drug disease (DrugCentral) file")
        i2_df = pd.read_csv(self.dc_ifile, sep="\t", header=0)
        i2_df = i2_df.replace(np.nan, '')  # replace NaN with empty string

        pubchemid = set(pubchemid)
        for i in range(i2_df.shape[0]):
            p = i2_df.iloc[i]['pubchem_cid']
            if p in pubchemid:
                continue
            else:
                pid = 'PUBCHEM_CID ' + str(p)
                cnodes['node_id'].append(pid)
                if "|" in dbid_dict[p]['names']:
                    cnodes['node_label'].append(dbid_dict[p]['names'].split("|")[0])
                else:
                    cnodes['node_label'].append(dbid_dict[p]['names'])
                cnodes['node_definition'].append(i2_df.iloc[i]['smiles'])
                cnodes['node_namespace'].append('IDG')
                cnodes['node_synonyms'].append(dbid_dict[p]['names'])
                cnodes['node_dbxrefs'].append(dbid_dict[p]['DrugBank'])
                cnodes['value'].append('')
                cnodes['lowerbound'].append('')
                cnodes['upperbound'].append('')
                cnodes['units'].append('')

        # write to an output dataframe
        print('write data to an output TSV file')
        o_df = pd.DataFrame(cnodes)
        o_df = o_df.drop_duplicates()
        o_df.to_csv(self.compound_node_file, sep="\t", index=False)

    def populate_disease_nodes(self):
        """
        Populate node fields and create a TSV file
        """
        print("\ngenerate disease nodes file")
        pnodes = {'node_id': [], 'node_namespace': [], 'node_label': [], 'node_definition': [], 'node_synonyms': [],
                  'node_dbxrefs': [], 'value': [], 'lowerbound': [], 'upperbound': [], 'units': []}

        if self.dc_ifile is None or self.disease_node_file is None:
            print("input/output file is missing")
            exit(-1)

        # read input file as dataframe and write data to output file
        print("fetch data from drug disease (DrugCentral) file")
        i_df = pd.read_csv(self.dc_ifile, sep="\t", header=0)
        i_df = i_df.replace(np.nan, '')  # replace NaN with empty string

        for col in i_df.columns:
            if col == 'snomed_conceptid':
                snomed_id = list(i_df['snomed_conceptid'])
                pnodes['node_id'] = ['SNOMEDCT_US ' + str(v) for v in snomed_id]
            elif col == 'omop_concept_name':
                pnodes['node_label'] = list(i_df['omop_concept_name'])
            elif col == 'omop_concept_id':
                omop_concept = list(i_df['omop_concept_id'])
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
        print('write data to an output TSV file')
        o_df = pd.DataFrame(pnodes)
        o_df = o_df.drop_duplicates()
        o_df.to_csv(self.disease_node_file, sep="\t", index=False)

    def populate_compound_disease_edges(self):
        """
        populate compound_disease edge fields and create a TSV file
        """
        print("\ngenerate compound disease edges file")
        # dict to store values
        edges = {'subject_id': [], 'relationship': [], 'object_id': [], 'evidence_class': []}

        if self.dc_ifile is None or self.compound_disease_edge_file is None:
            print("input/output file is missing")
            exit(-1)

        # read input file as dataframe and write data to output file
        print("fetch data from drug disease (DrugCentral) file")
        i_df = pd.read_csv(self.dc_ifile, sep="\t", header=0)
        i_df = i_df.replace(np.nan, '')  # replace NaN with empty string

        for col in i_df.columns:
            if col == 'snomed_conceptid':
                snomed_id = list(i_df['snomed_conceptid'])
                edges['object_id'] = ['SNOMEDCT_US ' + str(v) for v in snomed_id]
            elif col == 'pubchem_cid':
                pubchemid = i_df['pubchem_cid'].to_numpy(dtype=int)
                edges['subject_id'] = ['PUBCHEM_CID ' + str(v) for v in pubchemid]
            else:
                continue

        edges['relationship'] = ['indication'] * i_df.shape[0]
        edges['evidence_class'] = ['' for i in range(i_df.shape[0])]

        # write to an output dataframe
        print('write data to an output TSV file')
        o_df = pd.DataFrame(edges)
        o_df = o_df.drop_duplicates()
        o_df.to_csv(self.compound_disease_edge_file, sep="\t", index=False)

    def populate_compound_protein_edges(self):
        """
        populate compound_protein edge fields and create a TSV file
        """
        print("\ngenerate compound protein edges file")
        # dict to store values
        edges = {'subject_id': [], 'relationship': [], 'object_id': [], 'evidence_class': []}

        if self.tcrd_ifile is None or self.compound_protein_edge_file is None:
            print("input/output file is missing")
            exit(-1)

        # read input file as dataframe and write data to output file
        print("fetch data from compound activity (TCRD) file")
        i_df = pd.read_csv(self.tcrd_ifile, sep="\t", header=0)
        i_df = i_df.replace(np.nan, '')  # replace NaN with empty string

        for col in i_df.columns:
            if col == 'uniprot':
                uniprot = list(i_df['uniprot'])
                edges['object_id'] = ['UNIPROTKB ' + str(v) for v in uniprot]
            elif col == 'pubchem_cid':
                pubchemid = i_df['pubchem_cid'].to_numpy(dtype=int)
                edges['subject_id'] = ['PUBCHEM_CID ' + str(v) for v in pubchemid]
            elif col == 'act_type':
                edges['evidence_class'] = list(i_df['act_type'])
            else:
                continue

        edges['relationship'] = ['bioactivity'] * i_df.shape[0]

        # write to an output dataframe
        print('write data to an output TSV file')
        o_df = pd.DataFrame(edges)
        o_df = o_df.drop_duplicates()
        o_df.to_csv(self.compound_protein_edge_file, sep="\t", index=False)
