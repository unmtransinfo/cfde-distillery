// nodes - load data from compound file
LOAD CSV WITH HEADERS FROM "file:///cnodes.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Compound {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

// create index for compound nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (c:Compound) on (c.node_id)


// nodes - load data from protein file
LOAD CSV WITH HEADERS FROM "file:///pnodes.tsv" AS row FIELDTERMINATOR '\t'
CREATE (p:Protein {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

// create index for protein nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (p:Protein) on (p.node_id)


// edges (relationship) - load edges
:auto USING PERIODIC COMMIT 500
LOAD CSV WITH HEADERS FROM "file:///edges.tsv" AS row FIELDTERMINATOR '\t'
MATCH (c:Compound {node_id: trim(row.subject_id)}), (p:Protein {node_id: trim(row.object_id)})
CREATE (c)-[:RELATED_TO {relationship: row.relationship, evidence_class: row.evidence_class}]-> (p)
