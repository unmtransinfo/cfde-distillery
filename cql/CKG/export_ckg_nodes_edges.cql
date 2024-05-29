/*
After the export of nodes and edges, the respective files need 
to be copied to a directory and linked to docker container via
volume

Currently:

docker run --rm --volume /home/www/htdocs/x/cfde/distillery/data/CKG/data:/data \ 
--volume=/home/www/htdocs/x/cfde/distillery/data/CKG/import:/var/lib/neo4j/import \
--volume=/home/www/htdocs/x/cfde/distillery/data/CKG/conf:/var/lib/neo4j/conf \
--publish=3000:7474 --publish=7687:7687 \
--env NEO4J_PLUGINS='["apoc"]' \
-d neo4j:5.16.0

In this case, 
/home/www/htdocs/x/cfde/distillery/data/CKG/data, and
/home/www/htdocs/x/cfde/distillery/data/CKG/import are created, and the nodes
and edges file is placed inside the import folder
*/


/* Query to capture all nodes according to the SAB required
 and export as file (nodes.tsv). Output also contains a 
 column called Node_label, which stores the node label 
 (compound, disease, tissue etc.)
 Takes <2min to execute */

WITH ['EFO', 'GO', 'LINCS', 'NCBI', 'SNOMEDCT_US', 'PUBCHEM', 'GTEXEXP', 'EXPBINS', 'IDGD', 'UBERON', 'IDGP', 'UNIPROTKB', 'HGNC'] as SAB_LIST, 
"MATCH (t:Term)-[r1]-(cde:Code)-[r2]-(con:Concept) 
    WHERE cde.SAB in $sabList
    RETURN con.CUI as CUI, 
    CASE cde.SAB
        WHEN 'EFO' THEN 'EFO'
        when 'NCBI' then 'NCBI'
        when 'UNIPROTKB' then 'Protein'
        when 'HGNC' then 'Gene'
        when 'GO' then 'GO'
        when 'EFO' then 'EFO'
        when 'SNOMEDCT_US' then 'Disease'
        when 'PUBCHEM' then 'Compound'
        when 'EXPBINS' then 'Tissue'
        when 'UBERON' then 'Tissue'
        when 'GTEXEXP' then 'Tissue'
    end as Node_label,
    cde.SAB as SAB, cde.CODE as Code, t.name as Label, t.name as Synonyms" AS query
CALL apoc.export.csv.query(query,"nodes.tsv", {params: {sabList: SAB_LIST}, separator: '\t'}) 
YIELD file, nodes, relationships, properties, time
RETURN file, nodes, relationships, properties, time


/* Query to capture all relationships where the source and the
 target nodes' SAB belongs the required list
 Takes <20 min to execute (for batch size 100000), <18 min (for batch size 1000000) */

WITH ['EFO', 'GO', 'LINCS', 'NCBI', 'SNOMEDCT_US', 'PUBCHEM', 'GTEXEXP', 'EXPBINS', 'IDGD', 'UBERON', 'IDGP', 'UNIPROTKB', 'HGNC'] as SAB_LIST, 
"MATCH (source_code:Code)-[r1]-(source_concept:Concept)-[relationOfInterest]-(target_concept:Concept)-[r2]-(target_code:Code)
WHERE source_code.SAB IN $sabList AND target_code.SAB IN $sabList
RETURN source_concept.CUI as source, type(relationOfInterest) as relation, target_concept.CUI as target, source_code.CodeID as source_label, target_code.CodeID as target_label, relationOfInterest.SAB as SAB" as query
CALL apoc.export.csv.query(query, "relationships.tsv", {params: {sabList: SAB_LIST, batchSize: 100000, separator: '\t'}})
YIELD file, nodes, relationships, properties, time
RETURN file, nodes, relationships, properties, time
