import os 
import pandas as pd

def import_user_config():
    base_dir = os.path.abspath(os.path.dirname(__file__))
    config_file = os.path.join(base_dir, "user_args.txt")
    variables = {}

    with open(config_file) as f:
        for line in f:
            line = line.strip()
            if not line or '=' not in line:
                continue
            # line = line.replace("export ", "", 1)
            key, value = map(str.strip, line.split("=", 1))
            variables[key] = value.strip('"')

    return variables



def import_softwares():
    base_dir = os.path.abspath(os.path.dirname(__file__))
    config_file = os.path.join(base_dir, "software_location.txt")
    variables = {}

    with open(config_file) as f:
        for line in f:
            line = line.strip()
            if not line or '=' not in line:
                continue
            # line = line.replace("export ", "", 1)
            key, value = map(str.strip, line.split("=", 1))
            variables[key] = value.strip('"')

    return variables


def read_exclude_sampleids():
    base_dir = os.path.abspath(os.path.dirname(__file__))  # Get the absolute path of the current script
    exclude_sampleids_filepath = os.path.join(base_dir, 'exclude_sampleids.csv')
    exclude_sampleids_df = pd.read_csv(exclude_sampleids_filepath, dtype={"Sample_ID":str})
    exclude_sampleids_list = exclude_sampleids_df.Sample_ID.tolist()

    return exclude_sampleids_list

