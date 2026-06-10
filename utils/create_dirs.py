import os
import sys
import pandas as pd
from utils.validate_paths import validate_paths  # Custom function to check if provided paths exist
from utils.print_aesthetics import *

def create_round_dirs(round_dir, code, structure_type="output"):
    """
    Creates a directory structure for a specific 'round' (analysis run).
    
    Parameters:
    - round_dir (str): The base directory for the round.
    - dir_structure (str): Path to a text file listing subdirectories to create.
    - structure_type (str): Type of structure to create under the round_dir (either 'logs' or 'output').
    """

    # Ensure both round_dir and dir_structure path exist
    if structure_type == 'logs':
        dir_structure = f"{code}/user/round_logdir_structure.txt"
    elif structure_type == 'output':
        dir_structure = f"{code}/user/round_outdir_structure.txt"        
    else: 
        print(f"{FAIL} ERROR: Invalid structure_type '{structure_type}'. Must be 'logs' or 'output'.")
        sys.exit(1)
    validate_paths([round_dir, dir_structure])

    # Read the subdirectory names from the dir_structure file
    with open(dir_structure, 'r') as file:
        dirlist = [line.strip() for line in file.readlines()]

    # Create each subdirectory under the specified structure type
    inaccessible_directories = False

    for dirname in dirlist:
        dirpath = f"{round_dir}/{structure_type}/{dirname}"
        try:
            os.makedirs(dirpath, exist_ok=True)
            print(f"{PASS} PASS: Created {dirpath}")             
        except Exception as e:
            print(f"{FAIL} ERROR: Could not create {dirpath}")
            inaccessible_directories = True

    if inaccessible_directories:
        sys.exit(1)


def create_flowcell_dirs(metadata_filepath, flowcells_directory, genome_shorthand, code, log_dir, current_round):
    """
    Creates a nested directory structure for each unique flowcell ID found in metadata.

    Parameters:
    - metadata_filepath (str): CSV file containing flowcell metadata (must include 'flowcell_id' column).
    - flowcells_directory (str): Base directory where flowcell directories will be created.
    - genome_shorthand (str): Subfolder under each flowcell ID, typically representing a genome name or reference.
    - dir_structure (str): Path to a text file listing subdirectories to create under each flowcell/genome path.
    """

    # Validate input file and directory paths
    validate_paths([metadata_filepath, flowcells_directory, code, log_dir])

    # Extract unique flowcell IDs from the metadata
    df = pd.read_csv(metadata_filepath, dtype={'flowcell_id': str})
    unique_flowcell_ids = df.flowcell_id.unique().tolist()
    # print(f"unique_flowcell_ids: {sorted(unique_flowcell_ids)}")
    
    # save flowcell IDs to log file
    with open(os.path.join(log_dir, f'{current_round}_unique_flowcell_ids.csv'), 'w') as f:
        for fc in sorted(unique_flowcell_ids):
            f.write('%s\n' % fc)
    
    # Read the subdirectory names from the file
    dir_structure = f"{code}/user/flowcell_outdir_structure.txt"      
    with open(dir_structure, 'r') as file:
        dirlist = [line.strip() for line in file.readlines()]

    # Create directories for each flowcell ID
    inaccessible_directories = False

    for flowcell_id in unique_flowcell_ids:

        # Try creating topmost directory for flowcell ID
        dirpath = f"{flowcells_directory}/{flowcell_id}/{genome_shorthand}"
        try:
            os.makedirs(dirpath, exist_ok=True)
            print(f"{PASS} PASS: For flowcell_id: {flowcell_id}, created {dirpath}")            
        except Exception as e:
            print(f"{FAIL} ERROR: For flowcell_id: {flowcell_id}, could not create {dirpath}")
            inaccessible_directories = True            

        # Try creating subdirectories for flowcell ID
        for dirname in dirlist:
            dirpath = f"{flowcells_directory}/{flowcell_id}/{genome_shorthand}/{dirname}"
            try:
                os.makedirs(dirpath, exist_ok=True)
                # print(f"{PASS} PASS: For flowcell_id: {flowcell_id}, created {dirpath}")
            except Exception as e:
                print(f"{FAIL} ERROR: For flowcell_id: {flowcell_id}, could not create {dirpath}")
                inaccessible_directories = True

    if inaccessible_directories:
        sys.exit(1)
