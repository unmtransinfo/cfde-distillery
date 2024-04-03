MATCH (d:Disease)-[]-(c:Compound)-[]-(g1:Gene) WHERE d.CUI='C0917798'
MATCH (g1:Gene)-[]-(e:EFO)-[]-(g2:Gene)
RETURN DISTINCT g1.CUI as positive_genes_cui, g1.node_label as positive_genes_label, e.CUI as features, c.CUI as features_to_remove, g2.CUI as unknown_genes_cui, g2.node_label as unknown_genes_label
UNION
MATCH (d:Disease)-[]-(c:Compound)-[]-(g1:Gene) WHERE d.CUI='C0917798'
MATCH (g1:Gene)-[]-(e:Compound)-[]-(g2:Gene)
RETURN DISTINCT g1.CUI as positive_genes_cui, g1.node_label as positive_genes_label, e.CUI as features, c.CUI as features_to_remove, g2.CUI as unknown_genes_cui, g2.node_label as unknown_genes_label

