### <ins>Illuminating the Druggable Genome (IDG)</ins>

### Overview of the Illuminating the Druggable Genome (IDG) Tutorial

## Introduction

This tutorial provides a comprehensive guide to exploring the Illuminating the Druggable Genome (IDG) project using the Unified Biological Knowledge Graph (UBKG). The IDG program aims to enhance our understanding of the under-explored proteins in the druggable genome. By leveraging the vast data repository of UBKG, this tutorial demonstrates a series of Cypher query examples that allow users to unravel complex relationships between genes, diseases, and drugs.

### What You Will Learn

Through this tutorial, you will gain insights into:
<ul>
<li />Executing advanced queries in UBKG to map diseases and drugs.

<li />Analyzing gene expression profiles and their association with diseases.

<li />Exploring bioactivity relationships between compounds and proteins.

<li />Integrating multi-dimensional data from diverse datasets like LINCS, GTEx, and DrugCentral.
</ul>
 
### Target Audience

The content is primarily intended for bioinformaticians, researchers in genomics and pharmacology, and anyone interested in the intersection of big data and drug discovery.

### Tutorial Structure

The tutorial is structured into three main sections, each containing step-by-step examples:
<ol>
<li />**Introduction to IDG Queries within UBKG**: This section showcases how to map disease-related terms and compounds, using PUBCHEM and SNOMEDUS_CT as examples.
<li />**IDG Use-Case Combining IDG, LINCS, and GTEx Datasets**: Here, you'll learn to identify gene expression patterns and the impact of compounds on specific genes, integrating data from various sources.
<li />**Additional Use-Case Development**: Focused on the human ALOX genes and their role in asthma, this part dives into the specifics of querying genes and compounds related to inflammation and asthma.
</ol>
 
### Prerequisites

To make the most of this tutorial, familiarity with Neo4j, Cypher Query Language, and a basic understanding of genomics is recommended. Each example is accompanied by a detailed query and visual results to aid comprehension.

### Conclusion

By the end of this tutorial, you will be equipped with the knowledge to effectively utilize UBKG for exploring the intricate network of genes, diseases, and drugs, contributing significantly to the fields of genomics and drug discovery.
<br /><br />
Let's get started!

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

<strong>Example 3: Additional use case development inspired by the human ALOX genes encoding lipoxygenases and the inflammation associated with asthma</strong><br /><br />
<strong>Example 3a:</strong> Combining data from IDG-DrugCentral and LINCS, below are the results for a query seeking genes and Pubchem compounds associated with the disease “asthma”

```cypher 
WITH ['Asthma'] as theDisease
MATCH (hpo_concept:Concept)-[:CODE]-(hpo_code:Code {SAB:'HPO'})-[:PT]-(hpo_term:Term)
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
MATCH (pubchem_concept:Concept)-[:CODE]->(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH gr=(hpo_concept)-[r1:associated_with]-(gene_concept)-[r2 {SAB:'LINCS'}]-(pubchem_concept)
WHERE hpo_term.name in theDisease
RETURN *
LIMIT 50
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/3a.png?raw=true" width="100%">

A central node in the graph resulting from the above query, “C1412361,” is a gene named ALOX5 which is known to be involved with inflammation. Inflammation can cause swelling of the inner lining of the airways, leading to asthma symptoms. Specifically, human ALOX5 (Arachidonate 5-Lipoxygenase) is a member of the lipoxygenase (LOX) gene family and encodes an enzyme that performs the first two steps of leukotriene synthesis from arachidonic acid. ALOX5 is expressed in multiple cell types including bone marrow-derived cells, epithelial cells, skin cells, and cells involved in the regulation of either immune responses or inflammation. The leukotrienes synthesized by ALOX5 are known mediators of several allergic and inflammatory diseases, including asthma. The resulting graph shows that ALOX5 is negatively regulated by 39 compounds and positively regulated by 4 compounds within UBKG.


<strong>Example 3b: The human genome has six functional lipoxygenase genes (ALOX15, ALOX15B, ALOX12, ALOX12B, ALOXE3, and ALOX5). The query developed below searches for isoforms of lipoxygenase enzymes within UBKG containing the string “ALOX” and returns a sub-graph of enzymes and compounds associated via bioactivity:</strong><<br /><br />

```cypher 
MATCH (c_term:Term)-[:PREF_TERM]-(c_concept:Concept)-[:bioactivity {SAB: 'IDGP'}]-(p_concept:Concept)-[:PREF_TERM]-(p_term:Term)
WHERE p_term.name CONTAINS "ALOX"
RETURN c_concept, p_concept
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/3b.png?raw=true" width="100%">

This query returns several of the six known human lipoxygenase genes along with Pubchem compounds associated via bioactivity. The resulting graph shows there are several compounds associated with two different ALOX isoforms. Specifically, ten compounds are associated with both ALOX15 and ALOX12, and one compound is associated with both ALOX15 and ALOX15B.

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/3b_table.png?raw=true" width="75%">

Taking a closer look at one of the compounds associated with ALOX15B via bioactivity (PUBCHEM:1778842). The compound name is 3-[(4-Methylphenyl)methylsulfanyl]-1-phenyl-1,2,4-triazole. This compound, also known as “MLS000536924” was identified by Jameson II et al (2014) using high throughput screening as a competitive inhibitor (ki = 2.5+/-0.5 uM) of human epithelial 15-lipoxygenase-2 (ALOX15B). Several other compounds from the knowledge graph were manually investigated and found to be associated with Lipoxygenase activity in literature searches, supporting the results of our Cypher queries on UBKG.


<strong>Example 3c:</strong> In order to investigate the unique tissue types associated with the HGNC term “ALOX5”, we devise a query that combines data from both the IDG and the GTEx DCCs:

```cypher 
MATCH (uniprot_cui)-[:gene_product_of]->(hgnc_cui:Concept)-[:CODE]->(hgnc_code:Code {SAB:'HGNC'})-[]->(hgnc_term:Term {name:'ALOX5'})
MATCH (hgnc_cui)-[:expresses]->(gtexexp_cui:Concept)-[:CODE]->(gtexexp_code:Code {SAB:'GTEXEXP'})
MATCH (expbins_code:Code {SAB:'EXPBINS'})<-[:CODE]-(expbins_cui:Concept)-[:has_expression]-(gtexexp_cui)-[:expressed_in]->(ub_cui:Concept)-[:CODE]->(ub_code:Code {SAB:'UBERON'})-[:PT]-(ub_term:Term)
MATCH (ub_cui)-[:PREF_TERM]->(ub_term:Term)
MATCH (pubchem_cui2)-[:PREF_TERM]-(pubchem_2_term:Term)
RETURN distinct ub_term.name
```

Neo4j screenshot of query results:

<img src="https://github.com/unmtransinfo/cfde-distillery/blob/main/doc/UserGuide/images/3c.png?raw=true" width="100%">

The above query returns several locations of human cells expressing ALOX5 including breast epithelium, Peyer's patch, and anterior lingual gland. These results are consistent with a literature search focused on human ALOX5 expression profiles.
























