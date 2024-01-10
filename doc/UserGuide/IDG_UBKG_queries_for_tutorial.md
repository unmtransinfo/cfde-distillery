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

<strong>Example 2: IDG use-case which combines our IDG dataset with both LINCS and GTEx</strong><br /><br />
<strong>Example 2a:</strong> Find the top 25% genes that are highly expressed in the GTEx dataset using HGNC for gene annotations.

```cypher 
MATCH p=(tissue_code:Code {SAB:"GTEXEXP"})<-[:CODE]-(tissue_concept:Concept)-[r:expressed_in {SAB:"GTEXEXP"}]-(gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
WITH gene_code.CodeID as genes, COUNT(tissue_code.CodeID) as tissue_count, toInteger(COUNT(p) * 0.25) as top25Percent
ORDER BY tissue_count DESC
RETURN genes, tissue_count LIMIT 10
```

Neo4j screenshot of query results:
<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/2a.png?raw=true" width="100%">

<strong>Example 2b:</strong> Continuing from example 2a (where we found all genes in UBKG that are highly expressed in the GTEx dataset), we next focus on those genes which may be perturbed by a specific PubChem compound. This multi-DCC query relies on data from both the LINCS L1000 dataset and known drug targets found in data curated by IDG-DrugCentral. The first query below shows the result when all genes are considered: 

```cypher 
MATCH (tissue_concept:Concept)-[:CODE]->(tissue_code:Code {SAB:"GTEXEXP"})
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code {SAB:HGNC'})
MATCH (pubchem_concept:Concept)-[:CODE]->(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH (protein_concept:Concept)-[:CODE]->(protein_code:Code {SAB:"UNIPROTKB"})
MATCH (tissue_concept)-[r1:expressed_in {SAB:"GTEXEXP"}]-(gene_concept)-[r2 {SAB:LINCS'}]-(pubchem_concept)-[r3:bioactivity {SAB:'IDGP'}]-(protein_concept)
RETURN * LIMIT 5
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/2b_1.png?raw=true" width="100%">

And considering only the top 25% of genes:

```cypher 
MATCH p=(tissue_code:Code {SAB:"GTEXEXP"})<-[:CODE]-(tissue_concept:Concept)-[r:expressed_in {SAB:"GTEXEXP"}]-(gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
WITH gene_code.CodeID as genes, COUNT(tissue_code.CodeID) as tissue_count, toInteger(COUNT(p) * 0.25) as top25Percent
ORDER BY tissue_count DESC
UNWIND genes[..top25Percent] AS gene_id
WITH COLLECT(DISTINCT gene_id) as selectedGenes
MATCH (tissue_concept:Concept)-[:CODE]->(tissue_code:Code {SAB:"GTEXEXP"})
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code {SAB:HGNC'})
MATCH (pubchem_concept:Concept)-[:CODE]->(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH (protein_concept:Concept)-[:CODE]->(protein_code:Code {SAB:"UNIPROTKB"})
MATCH gr=(tissue_concept)-[r2:expressed_in {SAB:"GTEXEXP"}]-(gene_concept)-[r3 {SAB:'LINCS'}]-(pubchem_concept)-[r4:bioactivity {SAB:'IDGP'}]-(protein_concept)
WHERE gene_code.CodeID IN selectedGenes
RETURN gr LIMIT 5
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/2b_2.png?raw=true" width="100%">

<strong>Example 2c:</strong> Considering data from IDG-DrugCentral and LINCS, we next find all genes that are perturbed by compounds present in the LINCS dataset where those compounds are also linked to Parkinson's Disease. We take advantage of the “indication” relationship from DrugCentral to link the IDG-DrugCentral and LINCS datasets Below is the corresponding query:

```cypher 
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
MATCH (snomed_concept:Concept)-[:CODE]-(snomed_code:Code {SAB:'SNOMEDCT_US'})-[:PT]-(snomed_term:Term)
MATCH (pubchem_concept:Concept)-[:CODE]-(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH (gene_concept)-[r1 {SAB:'LINCS'}]-(pubchem_concept)-[r2:indication {SAB:'IDGD'}]-(snomed_concept)
WHERE snomed_term.name="Parkinson's disease"
RETURN * LIMIT 5
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/2c.png?raw=true" width="100%">

<strong>Example 2d:</strong> Find genes and compounds associated with the birth defect “Congenital diaphragmatic hernia”

```cypher 
WITH ['Congenital diaphragmatic hernia'] as birthDefects
MATCH (hpo_concept:Concept)-[:CODE]-(hpo_code:Code {SAB:'HPO'})-[:PT]-(hpo_term:Term)
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
MATCH (pubchem_concept:Concept)-[:CODE]->(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH gr=(hpo_concept)-[r1:associated_with]-(gene_concept)-[r2 {SAB:'LINCS'}]-(pubchem_concept)
WHERE hpo_term.name in birthDefects
RETURN *
LIMIT 10
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/2d.png?raw=true" width="100%">






