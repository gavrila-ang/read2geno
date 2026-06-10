#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 3 Check: Quality control of metadata file
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Cleaned metadata file
    - Sample processing input file
    - Raw, demultiplexed sample-level FASTQs

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Sample processing array file

Assumptions:
    - Custom input file directory structure

"""

##############################################################################################################################################################################

import os
import user
import utils
from time import time
from utils.print_aesthetics import *

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Declare input/output filepaths ===
clean_metadata_filepath = f"{round_directory}/output/01_flowcell_setup/flowcell_sheets/{current_round}_metadata_cleaned.csv"
sample_processing_input_filepath=f"{round_directory}/output/03_sample_processing/{current_round}_sample_processing_input.csv"
sample_processing_array_filepath=f"{round_directory}/output/03_sample_processing/03_sample_processing_arrays.csv"

# === Declare data type of special columns ===
dtype_mapping={'rfid':str}

##############################################################################################################################################################################

# === Runner ===

STEP = 0
elapsed_total = 0

def run_step(name, func):
    global STEP
    STEP += 1
    print(f"\n[{STEP:02d}] {name}")
    start = time()
    func()

    # Calculate time in hours and minutes
    elapsed = time() - start

def sh(cmd):
    """
    cmd: list[str]
    """
    subprocess.run(cmd, check=True)

##############################################################################################################################################################################

def main():

    print(f"Running sample processing checks....")

    # === Check if all the samples fastqs have been cleaned, aligned, and quality-controlled ===
    run_step(
        'Create another file containing arrays of sample_ids to be processed',
        lambda: utils.create_rerun_arrays.create_sample_processing_arrays(
            metadata_filepath=clean_metadata_filepath,
            sample_processing_input_filepath=sample_processing_input_filepath,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            exclude_sampleids_list=exclude_sampleids_list,
            output_filepath=sample_processing_array_filepath,
            dtype_mapping=dtype_mapping
            )
    )

    # === If the sample processing array filepath IS NOT empty, rerun Step 3 ===
    with open(sample_processing_array_filepath) as fp:
        lines = sum(1 for line in fp)

    if os.path.exists(sample_processing_array_filepath) and lines > 1:
        print("{FAIL} QC FAIL! There are samples that have not been successfully processed.")

    # === If the sample processing array filepath IS empty, create the input files for Step 3 ===
    else:
        print(f"{PASS} QC PASS! All sample_ids have been successfully processed.")

##############################################################################################################################################################################

if __name__ == "__main__":
    main()