// find node with Alzheimer's disease
MATCH (n:Disease) WHERE n.node_label CONTAINS 'Alzheimer' RETURN n;

// Find all compounds related to Alzheimer's disease
MATCH (d:Disease{node_label:"Alzheimer's disease"})-[r:RELATED_TO]-(c:Compound) RETURN c.node_id, c.node_label, c.node_dbxrefs;

// Find all compounds, proteins related to Alzheimer's disease. return graph
MATCH (d:Disease{node_label:"Alzheimer's disease"})-[r:RELATED_TO]-(c:Compound)-[r2:RELATED_TO]-(p:Protein) RETURN d,c,p;

// Find all compounds, proteins related to Alzheimer's disease. return some fields
MATCH (d:Disease{node_label:"Alzheimer's disease"})-[r:RELATED_TO]-(c:Compound)-[r2:RELATED_TO]-(p:Protein) RETURN p.node_id, p.node_synonyms, p.node_dbxrefs;
