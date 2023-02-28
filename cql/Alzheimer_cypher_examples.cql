MATCH (n:Disease) WHERE n.node_label CONTAINS 'Alzheimer' RETURN n;

MATCH (d:Disease{node_label:"Alzheimer's disease"})-[r:RELATED_TO]-(c:Compound) RETURN c.node_id, c.node_label, c.node_dbxrefs;

MATCH (d:Disease{node_label:"Alzheimer's disease"})-[r:RELATED_TO]-(c:Compound)-[r2:RELATED_TO]-(p:Protein) RETURN d,c,p;

MATCH (d:Disease{node_label:"Alzheimer's disease"})-[r:RELATED_TO]-(c:Compound)-[r2:RELATED_TO]-(p:Protein) RETURN p.node_id, p.node_synonyms, p.node_dbxrefs;
