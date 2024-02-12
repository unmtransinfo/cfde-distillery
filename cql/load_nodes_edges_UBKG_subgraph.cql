### All CSV files can be downloaded from the Chiltepin server (https://chiltepin.health.unm.edu/x/cfde/distillery/data/UBKG_subgraph.zip). After creating the Neo4j database, copy all CSV files to the import folder of the database. ###
  
// nodes - load data from compound file
LOAD CSV WITH HEADERS FROM "file:///PUBCHEM.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Compound {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids. It helps in loading edges quickly
MATCH (c:Compound)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index for compound nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (c:Compound) on (c.CUI);

------------------------------------

// nodes - load data from disease file
LOAD CSV WITH HEADERS FROM "file:///SNOMEDCT_US.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Disease {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids
MATCH (c:Disease)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index 
CREATE INDEX for (c:Disease) on (c.CUI);

------------------------------------

// nodes - load data from gene file
LOAD CSV WITH HEADERS FROM "file:///HGNC.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Gene {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids
MATCH (c:Gene)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index 
CREATE INDEX for (c:Gene) on (c.CUI);

------------------------------------

// nodes - load data from the Protein file
LOAD CSV WITH HEADERS FROM "file:///UNIPROTKB.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Protein {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids
MATCH (c:Protein)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index 
CREATE INDEX for (c:Protein) on (c.CUI);

------------------------------------

// nodes - load data from tissue file
LOAD CSV WITH HEADERS FROM "file:///EXPBINS.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Tissue {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

LOAD CSV WITH HEADERS FROM "file:///UBERON.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Tissue {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

LOAD CSV WITH HEADERS FROM "file:///GTEXEXP.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:Tissue {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids
MATCH (c:Tissue)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index 
CREATE INDEX for (c:Tissue) on (c.CUI);

------------------------------------

LOAD CSV WITH HEADERS FROM "file:///EFO.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:EFO {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids
MATCH (c:EFO)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index 
CREATE INDEX for (c:EFO) on (c.CUI);

------------------------------------

LOAD CSV WITH HEADERS FROM "file:///GO.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:GO {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids
MATCH (c:GO)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index 
CREATE INDEX for (c:GO) on (c.CUI);

------------------------------------

LOAD CSV WITH HEADERS FROM "file:///NCBI.tsv" AS row FIELDTERMINATOR '\t'
CREATE (c:NCBI {CUI:trim(row.CUI), SAB:trim(row.SAB), node_code:trim(row.Code), node_label:trim(row.Label), value:trim(row.Value), lower_bound:trim(row.Lower_bound), upper_bound:trim(row.Upper_bound), unit:trim(row.Unit), node_synonyms:trim(row.Synonyms)});

// delete duplicate node ids. It helps in loading edges quickly
MATCH (c:NCBI)
WITH c.CUI AS cui, COLLECT(c) AS cnodes
WHERE SIZE(cnodes) > 1
FOREACH (n IN TAIL(cnodes) | DETACH DELETE n);

// create index for compound nodes on node_id. It helps in loading edges quickly
CREATE INDEX for (c:NCBI) on (c.CUI);

------------------------------------
// Run the following CQL for all possible combinations of node labels. The following CQL uses Compound and Disease.
:auto USING PERIODIC COMMIT 500
LOAD CSV WITH HEADERS FROM "file:///edges_pk.tsv" AS row FIELDTERMINATOR '\t'
MATCH (c1: Compound {CUI: trim(row.source)}), (c2: Disease {CUI: trim(row.target)})
CALL apoc.create.relationship(c1, row.relation, {SAB: row.SAB, evidence_class: row.evidence}, c2) YIELD rel
return count(rel);
