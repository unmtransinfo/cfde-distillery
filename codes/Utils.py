import yaml

def parse_yml(file_path):
    with open(file_path) as r:
        try:
            valuesYaml = yaml.safe_load(r)
        except Exception as e:
            valuesYaml = None
        
    return valuesYaml