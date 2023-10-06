
<h1><strong>Exhibit 1: IDG’s assigned use case (combining our IDG dataset with LINCS and GTEx)</strong></h1>
<br />
<h2><i>These queries were executed on the UBKG instance first loaded to Chiltepin on September 13, 2023</i></h2>

<br />

<h2>Exhibit 1a: Find the top 25% genes that are highly expressed in the GTEx <tissue> dataset using HGNC for gene annotations.</h2>

````
MATCH p=(tissue_code:Code {SAB:"GTEXEXP"})<-[:CODE]-(tissue_concept:Concept)-[r:expressed_in {SAB:"GTEXEXP"}]-(gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
WITH gene_code.CodeID as genes, COUNT(tissue_code.CodeID) as tissue_count, toInteger(COUNT(p) * 0.25) as top25Percent
ORDER BY tissue_count DESC
RETURN genes, tissue_count LIMIT 10
````

<h2>Exhibit 1b: Continuing from exhibit 1a (where we found all genes in UBKG that are highly expressed in the GTEx <tissue> dataset), we next focus on those genes which may be perturbed by a specific PubChem compound.  This multi-DCC query relies on data from  both the LINCS L1000 dataset and known drug targets found in data curated by IDG.
The first query below shows the result when all genes are considered:
</h2>

````
MATCH (tissue_concept:Concept)-[:CODE]->(tissue_code:Code {SAB:"GTEXEXP"})
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code {SAB:HGNC'})
MATCH (pubchem_concept:Concept)-[:CODE]->(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH (protein_concept:Concept)-[:CODE]->(protein_code:Code {SAB:"UNIPROTKB"})
MATCH (tissue_concept)-[r1:expressed_in {SAB:"GTEXEXP"}]-(gene_concept)-[r2 {SAB:LINCS'}]-(pubchem_concept)-[r3:bioactivity {SAB:'IDGP'}]-(protein_concept)
RETURN * LIMIT 5
````

<h2>And considering only the top 25% of genes:</h2>


````
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
````

<h2>Exhibit 1c: Considering data from IDG and LINCS, we next find all genes that are perturbed by compounds present in the LINCS dataset where those compounds are also linked to Parkinson's Disease.  We take advantage of the “indication” relationship from DrugCentral to link the IDG and LINCS datasets  Below is the corresponding query:</h2>

````
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
MATCH (snomed_concept:Concept)-[:CODE]-(snomed_code:Code {SAB:'SNOMEDCT_US'})-[:PT]-(snomed_term:Term)
MATCH (pubchem_concept:Concept)-[:CODE]-(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH (gene_concept)-[r1 {SAB:'LINCS'}]-(pubchem_concept)-[r2:indication {SAB:'IDGD'}]-(snomed_concept)
WHERE snomed_term.name="Parkinson's disease"
RETURN * LIMIT 5
````

<h2>Exhibit 1d: Find genes and compounds associated with the birth defect “Congenital diaphragmatic hernia”</h2>

````
WITH ['Congenital diaphragmatic hernia'] as birthDefects
MATCH (hpo_concept:Concept)-[:CODE]-(hpo_code:Code {SAB:'HPO'})-[:PT]-(hpo_term:Term)
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
MATCH (pubchem_concept:Concept)-[:CODE]->(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH gr=(hpo_concept)-[r1:associated_with]-(gene_concept)-[r2 {SAB:'LINCS'}]-(pubchem_concept)
WHERE hpo_term.name in birthDefects
RETURN *
LIMIT 10

````




