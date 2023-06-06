// nodes - load data from compound file
LOAD CSV WITH HEADERS FROM "file:///nodes_compound.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Compound {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

// nodes - load data from LINCS compound file
LOAD CSV WITH HEADERS FROM "file:///lincs_nodes_chem.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Compound {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

// create index for compound nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (c:Compound) on (c.node_id)

// delete duplicate node ids.It helps in loading edges quickly
MATCH (c:Compound)
WITH c.node_id AS node_id, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

---
// nodes - load data from protein file
LOAD CSV WITH HEADERS FROM "file:///nodes_protein.tsv" AS row FIELDTERMINATOR '\t'
CREATE (p:Protein {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

// create index for protein nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (p:Protein) on (p.node_id)

// delete duplicate node ids.It helps in loading edges quickly
MATCH (p:Protein)
WITH p.node_id AS node_id, COLLECT(p) AS pnodes
WHERE SIZE(pnodes) > 1
FOREACH (n IN TAIL(pnodes) | DETACH DELETE n);

----
// nodes - load data from disease file
LOAD CSV WITH HEADERS FROM "file:///nodes_disease.tsv" AS row FIELDTERMINATOR '\t'
CREATE (d:Disease {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

// create index for disease nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (d:Disease) on (d.node_id);

// delete duplicate node ids.It helps in loading edges quickly
MATCH (d:Disease)
WITH d.node_id AS node_id, COLLECT(d) AS dnodes
WHERE SIZE(dnodes) > 1
FOREACH (n IN TAIL(dnodes) | DETACH DELETE n);

----
// nodes - load data from LINCS gene file
LOAD CSV WITH HEADERS FROM "file:///lincs_nodes_gene_1.tsv" AS row FIELDTERMINATOR '\t'
CREATE (g:Gene {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

LOAD CSV WITH HEADERS FROM "file:///lincs_nodes_gene_2.tsv" AS row FIELDTERMINATOR '\t'
CREATE (g:Gene {node_id:trim(row.node_id), node_namespace:trim(row.node_namespace), node_label:trim(row.node_label), node_definition
:trim(row.node_definition), node_synonyms:trim(row.node_synonyms), node_dbxrefs:trim(row.node_dbxrefs), value:trim(row.value), lowerbound
:trim(row.lowerbound), upperbound:trim(row.upperbound), units:trim(row.units)});

// create index for Gene nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (g:Gene) on (g.node_id)

// delete duplicate node ids.It helps in loading edges quickly
MATCH (g:Gene)
WITH g.node_id AS node_id, COLLECT(g) AS gnodes
WHERE SIZE(gnodes) > 1
FOREACH (n IN TAIL(gnodes) | DETACH DELETE n);

----
// bioactivity edges (relationship) - load edges
:auto USING PERIODIC COMMIT 500
LOAD CSV WITH HEADERS FROM "file:///edges_bioactivity_compound_protein.tsv" AS row FIELDTERMINATOR '\t'
MATCH (c:Compound {node_id: trim(row.subject_id)}), (p:Protein {node_id: trim(row.object_id)})
CREATE (c)-[:Bioactivity {relationship: row.relationship, evidence_class: row.evidence_class}]-> (p)

// indication edges (relationship) - load edges
:auto USING PERIODIC COMMIT 500
LOAD CSV WITH HEADERS FROM "file:///edges_indication_compound_disease.tsv" AS row FIELDTERMINATOR '\t'
MATCH (c:Compound {node_id: trim(row.subject_id)}), (d:Disease {node_id: trim(row.object_id)})
CREATE (c)-[:Indication {relationship: row.relationship, evidence_class: row.evidence_class}]-> (d)

// Regulated_by edges (relationship) - load edges
:auto USING PERIODIC COMMIT 500
LOAD CSV WITH HEADERS FROM "file:///lincs_edges_gene_1.tsv" AS row FIELDTERMINATOR '\t'
MATCH (c:Compound {node_id: trim(row.object)}), (g:Gene {node_id: trim(row.subject)})
CREATE (c)-[:Regulated_by {relationship: row.relationship, evidence_class: row.evidence_class, uri: row.uri}]-> (g)

// Regulated_by edges (relationship) - load edges
:auto USING PERIODIC COMMIT 500
LOAD CSV WITH HEADERS FROM "file:///lincs_edges_gene_2.tsv" AS row FIELDTERMINATOR '\t'
MATCH (c:Compound {node_id: trim(row.object)}), (g:Gene {node_id: trim(row.subject)})
CREATE (c)-[:Regulated_by {relationship: row.relationship, evidence_class: row.evidence_class, uri: row.uri}]-> (g)

// Similar_to edges (relationship) - load edges
:auto USING PERIODIC COMMIT 500
LOAD CSV WITH HEADERS FROM "file:///lincs_edges_chem.tsv" AS row FIELDTERMINATOR '\t'
MATCH (c1:Compound {node_id: trim(row.subject)}), (c2:Compound {node_id: trim(row.object)})
CREATE (c1)-[:Similar_to {relationship: row.relationship, evidence_class: row.evidence_class, uri: row.uri}]-> (c2)

