
<strong><h1>NOTE: these example queries were each executed successfully on the UBKG instance first loaded on chiltepin September 13, 2023</h1></strong>

<h2>Showing the IDGD (IDG-Disease) mapping between PUBCHEM and SNOMEDUS_CT:</h2>

````
MATCH (pubchem_code:Code {SAB:'PUBCHEM'})-[:CODE]-(pubchem_cui:Concept)-[:indication {SAB:'IDGD'}]-(snomed_cui:Concept)-[:CODE]-(snomed_code:Code {SAB:"SNOMEDCT_US"})
RETURN * LIMIT 10
````

<h2>Showing the results for disease terms containing the string “diabetes” and linked via the IDG indication relationship:</h2>

````
MATCH (t1:Term)<-[:PREF_TERM]-(p1:Concept)-[:indication]-(c2:Concept)--(t2:Term) WHERE t1.name 
CONTAINS 'diabetes' 
RETURN * ;
````

<h2>Listing compounds related to IDG terms containing the string “diabetes”:</h2>

````
MATCH (c_term:Term)-[:PREF_TERM]-(c_concept:Concept)-[:indication {SAB: 'IDGD'}]-(d_concept:Concept)-[:PREF_TERM]-(d_term:Term)
WHERE d_term.name CONTAINS "diabetes"
RETURN toLower(c_term.name) as Name
ORDER BY Name
````

<h2>Showing compounds and proteins (using SAB: IDG-P) related by bioactivity where the protein name contains the string “Cytochrome P450”:</h2>

````
MATCH (c_term:Term)-[:PREF_TERM]-(c_concept:Concept)-[:bioactivity {SAB: 'IDGP'}]-(p_concept:Concept)-[:PREF_TERM]-(p_term:Term)
WHERE p_term.name CONTAINS "Cytochrome P450"
RETURN *
````
