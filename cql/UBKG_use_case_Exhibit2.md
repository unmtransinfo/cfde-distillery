
<h1><strong>Exhibit 2: Additional use case development inspired by the human ALOX genes encoding lipoxygenases and the inflammation associated with asthma</strong></h1>

<h2><i>These queries were executed on the UBKG instance first loaded to Chiltepin on September 13, 2023</i></h2>

<br />

<h2>Exhibit 2a: Combining data from IDG and LINCS, below are the results for a query seeking genes and Pubchem compounds associated with the disease “asthma”</h2>

````
WITH ['Asthma'] as theDisease
MATCH (hpo_concept:Concept)-[:CODE]-(hpo_code:Code {SAB:'HPO'})-[:PT]-(hpo_term:Term)
MATCH (gene_concept:Concept)-[:CODE]->(gene_code:Code{SAB:'HGNC'})
MATCH (pubchem_concept:Concept)-[:CODE]->(pubchem_code:Code {SAB:'PUBCHEM'})
MATCH gr=(hpo_concept)-[r1:associated_with]-(gene_concept)-[r2 {SAB:'LINCS'}]-(pubchem_concept)
WHERE hpo_term.name in theDisease
RETURN *
LIMIT 50
````

<h2>A central node in the graph resulting from the above query, “C1412361,” is a gene named ALOX5 which is known to be involved with inflammation.  Inflammation can cause swelling of the inner lining of the airways, leading to asthma symptoms.  Specifically, human ALOX5 (Arachidonate 5-Lipoxygenase) is a member of the lipoxygenase (LOX) gene family and encodes an enzyme that performs the first two steps of leukotriene synthesis from arachidonic acid.  ALOX5 is expressed in multiple cell types including bone marrow-derived cells, epithelial cells, skin cells, and cells involved in the regulation of either immune responses or inflammation.  The leukotrienes synthesized by ALOX5 are known mediators of several allergic and inflammatory diseases, including asthma.  The resulting graph shows that ALOX5 is negatively regulated by 39 compounds and positively regulated by 4 compounds within UBKG.</h2>




<h2>Exhibit 2b: The human genome has six functional lipoxygenase genes (ALOX15, ALOX15B, ALOX12, ALOX12B, ALOXE3, and ALOX5).  The query developed below searches for isoforms of lipoxygenase enzymes within UBKG containing the string “ALOX” and returns a sub-graph of enzymes and compounds associated via bioactivity:</h2>

````
MATCH (c_term:Term)-[:PREF_TERM]-(c_concept:Concept)-[:bioactivity {SAB: 'IDGP'}]-(p_concept:Concept)-[:PREF_TERM]-(p_term:Term)
WHERE p_term.name CONTAINS "ALOX"
RETURN c_concept, p_concept
````

<h2>This query returns several of the six known human lipoxygenase genes along with Pubchem compounds associated via bioactivity. The resulting graph shows there are several compounds associated with two different ALOX isoforms. Specifically, ten compounds are associated with both ALOX15 and ALOX12, and one compound is associated with both ALOX15 and ALOX15B.</h2>

<h2>Taking a closer look at one of the compounds associated with ALOX15B via bioactivity (PUBCHEM:1778842).  The compound name is 3-[(4-Methylphenyl)methylsulfanyl]-1-phenyl-1,2,4-triazole.  This compound, also known as “MLS000536924” was identified by Jameson II et al (2014) using high throughput screening as a competitive inhibitor (ki = 2.5+/-0.5 uM) of human epithelial 15-lipoxygenase-2 (ALOX15B).  Several other compounds from the knowledge graph were manually investigated and found to be associated with Lipoxygenase activity in literature searches, supporting the results of our Cypher queries on UBKG.</h2>


<h2>Exhibit 2c: In order to investigate the unique tissue types associated with the HGNC term “ALOX5”, we devise a query that combines data from both the IDG and the GTEx DCCs:</h2>

````
MATCH (uniprot_cui)-[:gene_product_of]->(hgnc_cui:Concept)-[:CODE]->(hgnc_code:Code {SAB:'HGNC'})-[]->(hgnc_term:Term {name:'ALOX5'})
MATCH (hgnc_cui)-[:expresses]->(gtexexp_cui:Concept)-[:CODE]->(gtexexp_code:Code {SAB:'GTEXEXP'})
MATCH (expbins_code:Code {SAB:'EXPBINS'})<-[:CODE]-(expbins_cui:Concept)-[:has_expression]-(gtexexp_cui)-[:expressed_in]->(ub_cui:Concept)-[:CODE]->(ub_code:Code {SAB:'UBERON'})-[:PT]-(ub_term:Term)
MATCH (ub_cui)-[:PREF_TERM]->(ub_term:Term)
MATCH (pubchem_cui2)-[:PREF_TERM]-(pubchem_2_term:Term)
RETURN distinct ub_term.name
````

<h2>The above query returns several locations of human cells expressing ALOX5 including breast epithelium, Peyer's patch, and anterior lingual gland.  These results are consistent with a literature search focused on human ALOX5 expression profiles.</h2>





