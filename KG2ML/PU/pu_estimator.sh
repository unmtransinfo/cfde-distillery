set -x
time python TiCE/tice_pu_driver.py -inpfile /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/schizophrenia.csv -gene_label_file /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/gene_label_data.csv 
time python TiCE/tice_pu_driver.py -inpfile /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/raynaud.csv -gene_label_file /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/gene_label_data.csv 
time python TiCE/tice_pu_driver.py -inpfile /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/pulmonary_edema.csv -gene_label_file /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/gene_label_data.csv 
time python TiCE/tice_pu_driver.py -inpfile /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/parkinson.csv -gene_label_file /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/gene_label_data.csv 
time python TiCE/tice_pu_driver.py -inpfile /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/malaria.csv -gene_label_file /mnt/bigHDD/work/CFDE_all/CFDE_data/KG2ML/gene_label_data.csv 
set +x

