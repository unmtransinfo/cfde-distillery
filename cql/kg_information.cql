// number of nodes in KG

MATCH (n) RETURN COUNT(n);

// number of Protein nodes

MATCH (p:Protein) RETURN COUNT(p)


// number of Compound nodes

MATCH (c:Compound) RETURN COUNT(c)


// number of relationships in KG

MATCH p=()-->() RETURN COUNT(p)
MATCH p=()-[r:RELATED_TO]->() RETURN COUNT(p)


// find all unique evidence_class in the relationship

MATCH p=()-[r:RELATED_TO]->() WHERE r.evidence_class IS NOT NULL RETURN DISTINCT(r.evidence_class)


// find count of all evidence_class in the relationship

MATCH p=()-[r:RELATED_TO]->() WHERE r.evidence_class IS NOT NULL RETURN r.evidence_class AS evidence_class, COUNT(*) AS count


// find the number of relationships for a particular protein node

MATCH (p:Protein)-[r:RELATED_TO]-(c:Compound) WHERE p.node_label='CTH' RETURN COUNT(r)


