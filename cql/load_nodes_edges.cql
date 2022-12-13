// nodes - load data from compound file

LOAD CSV WITH HEADERS FROM "file:///cnodes.tsv"
AS row FIELDTERMINATOR '\t' WITH row
CREATE (c:Compound {node_id:row.node_id, node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});



// nodes - load data from protein file

LOAD CSV WITH HEADERS FROM "file:///pnodes.tsv"
AS row FIELDTERMINATOR '\t' WITH row
CREATE (p:Protein {node_id:row.node_id, node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});


// edges (relationship) - load edges

:auto USING PERIODIC COMMIT 100
LOAD CSV WITH HEADERS FROM "file:///edges.tsv" 
AS row FIELDTERMINATOR '\t' WITH row
MATCH (c:Compound {node_id:trim(row.subject_id)}), (p:Protein {node_id:trim(row.object_id)})
CREATE (c)-[:RELATED_TO{relationship: row.relationship, evidence_class: row.evidence_class}]->(p) ;

