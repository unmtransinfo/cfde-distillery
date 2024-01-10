### <ins>Illuminating the Druggable Genome (IDG)</ins>

<strong>Example 1: Introduction to IDG queries within UBKG</strong><br /><br />
<strong>Example 1a:</strong> Showing the IDGD (IDG-Disease) mapping between PUBCHEM and SNOMEDUS_CT:

```cypher 
MATCH (pubchem_code:Code {SAB:'PUBCHEM'})-[:CODE]-(pubchem_cui:Concept)-[:indication {SAB:'IDGD'}]-(snomed_cui:Concept)-[:CODE]-(snomed_code:Code {SAB:"SNOMEDCT_US"})
RETURN * LIMIT 10
```
Neo4j screenshot of query results:
<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/1a.png?raw=true" width="100%">

<strong>Example 1b:</strong> Showing the results for disease terms containing the string “diabetes” and linked via the IDG-DrugCentral indication relationship:

```cypher 
MATCH (t1:Term)<-[:PREF_TERM]-(p1:Concept)-[:indication]-(c2:Concept)--(t2:Term) WHERE t1.name 
CONTAINS 'diabetes' 
RETURN * ;
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/1b.png?raw=true" width="100%">

<strong>Example 1c:</strong> Listing compounds related to IDG terms containing the string “diabetes”:

```cypher 
MATCH (c_term:Term)-[:PREF_TERM]-(c_concept:Concept)-[:indication {SAB: 'IDGD'}]-(d_concept:Concept)-[:PREF_TERM]-(d_term:Term)
WHERE d_term.name CONTAINS "diabetes"
RETURN toLower(c_term.name) as Name
ORDER BY Name
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/1c.png?raw=true" width="100%">

<strong>Example 1d:</strong> Showing compounds and proteins (using SAB: IDG-P) related by bioactivity where the protein name contains the string “Cytochrome P450”:

```cypher 
MATCH (c_term:Term)-[:PREF_TERM]-(c_concept:Concept)-[:bioactivity {SAB: 'IDGP'}]-(p_concept:Concept)-[:PREF_TERM]-(p_term:Term)
WHERE p_term.name CONTAINS "Cytochrome P450"
RETURN *
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/1d.png?raw=true" width="100%">

