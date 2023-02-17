# cfde-distillery
CFDE Data Distillery Partnership Project IDG contributions

# nodes and edges files
All input and output files are here: https://drive.google.com/drive/folders/1eUcYVayYHM90ESrQqpx8GAcEYOXUCr-i and https://chiltepin.health.unm.edu/x/cfde/distillery/data/

# How to load data from TSV files to Neo4J KG
- Download compound_activity.tsv, drug_disease.csv, and drugbank_id.pkl from https://drive.google.com/drive/folders/1eUcYVayYHM90ESrQqpx8GAcEYOXUCr-i 
- If TCRD/DrugCentral databases changed, run compound_activity.sql and drug_central_data.sql present in the "sql" folder to get input data from TCRD and DrugCentral databases.
- Update the input and output files path in drugcentral_nodes.py, drugcentral_edges.py, nodes.py, and edges.py and run these codes.
- Move all output TSV files for nodes and edges to the "import" folder of the Neo4J database you created.
- Run all Cypher queries present in the "cql" folder.

