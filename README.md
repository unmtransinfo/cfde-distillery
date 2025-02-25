# cfde-distillery <img align="right" src="doc/images/cfde_logo_above_idg_logo.png" height="240">

## Overview

The CFDE Data Distillery Partnership Project, led by the CFDE HuBMAP, SenNet, and Kids
First teams, from Pitt and Children's Hospital of Philadelphia (CHOP) has developed a
Data Distillery Knowledge Graph (DDKG), which distills Common Fund data with semantic
interoperability, to support integrative biomedical research questions and 
scientific use cases.

The IDG-DCC team brings to this project a strong record of research in KG design,
analytics and KG-based ML. This repo represents the contributions of the IDG team
at the University of New Mexico. The IDG team has focused on including IDG datasets
into the DDKG, and developing IDG scientific use cases combining datasets from IDG and multiple
other Common Fund programs. In addition, the IDG team has developed the _Condensed KG_,
derived from the DDKG, but subsetted and scoped for IDG use cases, and simplified for
interpretability and usability, 

## DDKG Nodes and Edges files

All input and output files are here: https://drive.google.com/drive/folders/1eUcYVayYHM90ESrQqpx8GAcEYOXUCr-i and https://chiltepin.health.unm.edu/x/cfde/distillery/data/

### How to load data from TSV files to Neo4J KG
- Download compound_activity_input.tsv, drug_disease_input.tsv, and drugbank_id.pkl from https://drive.google.com/drive/folders/1eUcYVayYHM90ESrQqpx8GAcEYOXUCr-i
- If TCRD/DrugCentral databases changed, run compound_activity.sql and drug_central_data.sql present in the "sql" folder to get input data from TCRD and DrugCentral databases.
- Update the input and output files path in drugcentral_nodes.py, drugcentral_edges.py, nodes.py, and edges.py and run these codes.
- Move all output TSV files for nodes and edges to the "import" folder of the Neo4J database you created.
- In Neo4J database settings, set dbms.memory.heap.max_size=8G. If it gives error, increase the value to a bigger number (e.g. 16gb).
- Run all Cypher queries present in the "cql" folder.

## Condensed-KG

Hosting the comprehensive DDKG requires significant computational resources, and for IDG 
applications, a small fraction of the data is required, and a lightweight subset is advantageous.
To address this challenge, we developed a localized, cloud-ready, containerized KG.
Drawing inspiration from the LINCS-DCC (Maayan Lab) project, we
created a "Condensed-KG" derived and subsetted from the latest DDKG. This Condensed-KG,
consisting of 8 distinct node labels and 1042 distinct relationship types, is designed to simplify
and enhance interpretability, specifically tailored to IDG use cases.

<img src="doc/images/DDKG_to_CKG_workflow.png" height="180">

## Links

 * [NIH-CFDE.org](https://www.nih-cfde.org/)
 * [Distillery Partnership Home Page](https://github.com/nih-cfde/data-distillery/)
 * [Unified Biomedical Knowledge Graph (UBKG)](https://ubkg.docs.xconsortia.org/)
 * [CFDE_DataDistillery GitHub Repository](https://github.com/TaylorResearchLab/CFDE_DataDistillery)
 * [Poster: Illuminating the Druggable Genome (IDG) scientific use cases powered by the CFDE Data Distillery biomedical knowledge graph, integrating multiple Common Fund datasets](https://doi.org/10.5281/zenodo.10895777), Common Fund Data Ecosystem (CFDE) Meeting, Mar 19-20, 2024. 
