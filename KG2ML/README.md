# How to run gene_pu_classifier.py
To run this Python code, you need to provide: 1) an input file containing the CQL output for a given condition, 2) a file containing genes and their corresponding labels, and 3) a filename to save the bar plot. CQL for selecting positive and unknown genes and features for Pakinson's disease in the cql folder. Here is an example.
```
python gene_pu_classifier.py -inpfile ../../cfde-distillery_data/KG2ML/cerebral_all.csv -gene_label_file inputData/gene_label_data.csv -pltfile inputData/newplt.png
```
