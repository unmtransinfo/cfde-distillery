
/*
There might be a need to change the FIELDTERMINATOR in case of loading
nodes and edges. If there is an issue with current FIELDTERMINATOR, use
"," as the FIELDTERMINATOR
*/



// Query to load nodes in the database
// Uses APOC for batch mode

CALL apoc.periodic.iterate("LOAD CSV WITH HEADERS FROM 'file:///nodes.tsv' AS row FIELDTERMINATOR ','
    RETURN row", 
"WITH row, trim(row.Node_label) AS label

CALL apoc.create.node([label], {
    CUI: trim(row.CUI),
    SAB: trim(row.SAB),
    node_code: trim(row.Code),
    node_label: trim(row.node_label),
    node_synonyms: trim(row.Synonyms)
}) YIELD node

RETURN node;",
{batchSize:500000, parallel:true}
)
YIELD batches, total
RETURN batches,total

// create index for nodes based on CUI
MATCH (c)
WHERE c.CUI IS NOT NULL
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// delete duplicate nodes
CALL apoc.meta.nodeTypeProperties() YIELD nodeType
WITH COLLECT(nodeType) AS nodeLabels
UNWIND nodeLabels AS label
CALL apoc.schema.assert({label: ["CUI"]}, {}) YIELD label AS indexLabel
RETURN indexLabel;

// Query to create relationships in the database
CALL apoc.periodic.iterate(
  "LOAD CSV WITH HEADERS FROM 'file:///relationships.tsv' AS row FIELDTERMINATOR ','
   RETURN row",
  "MATCH (sourceNode {CUI: row.source})
   MATCH (targetNode {CUI: row.target})
   CALL apoc.create.relationship(sourceNode, row.relation, {SAB: row.SAB}, targetNode)
   YIELD rel
   RETURN rel",
  {batchSize:500000, parallel:true, iterateList: true, retries:5}
)


