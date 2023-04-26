- fetch data from DrugCentral database for "concept" node and "RELATED_TO" relationship

SELECT DISTINCT
  omop.concept_id omop_concept_id,
  omop.concept_name omop_concept_name,
  omop.snomed_conceptid, 
  ids.identifier pubchem_cid,
  s.id drug_id,
  s.name drug_name,
  s.smiles 
FROM
  omop_relationship omop
JOIN
  structures s ON omop.struct_id = s.id
JOIN
    identifier ids ON ids.struct_id = s.id
WHERE
  ids.id_type = 'PUBCHEM_CID'
  AND omop.relationship_name = 'indication'
  AND omop.snomed_conceptid IS NOT null
