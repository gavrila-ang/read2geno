#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 2 Check: Quality control of metadata file
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Cleaned metadata file
    - Demux input file
    - Demux array file
    - Raw, demultiplexed sample-level FASTQs 

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Demux array file
    - Sample processing input file

Assumptions:
    - Custom input file directory structure

"""

##############################################################################################################################################################################

import user
import utils
import os
from time import time

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Declare input/output filepaths ===
clean_metadata_filepath = f"{round_directory}/output/01_flowcell_setup/flowcell_sheets/{current_round}_metadata_cleaned.csv"
demux_input_filepath=f"{round_directory}/output/02_demux/{current_round}_demux_input.csv"
demux_array_filepath=f"{round_directory}/output/02_demux/02_demux_arrays.csv"
undemuxed_samples_filepath=f"{round_directory}/output/02_demux/undemuxed_samples.csv"
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

    print(f"Running demutiplexing checks....")

    # === Check if all the libraries have been demultiplexed successfully ===
    run_step(
        'Create a file containing arrays of libraries to be demuxed',
        lambda: utils.create_rerun_arrays.create_demux_arrays(
            demux_input_filepath=demux_input_filepath,
            metadata_filepath=clean_metadata_filepath,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            exclude_sampleids_list=exclude_sampleids_list,
            summary_filepath=undemuxed_samples_filepath,
            output_filepath=demux_array_filepath,
            dtype_mapping=dtype_mapping
        )
    )

    # === If the demux array filepath IS NOT empty, rerun Step 2 ===
    if os.path.exists(demux_array_filepath) and os.path.getsize(demux_array_filepath) != 0:
        print("There are libraries that have not been successfully demultiplexed.")

    # === If the demux array filepath IS empty, create the input files for Step 3 ===
    else:
        print("All libraries have been successfully demultiplexed")
        print(f"Creating sample processing input and array files...")

        run_step(
            'Create a file containing sample_ids and their processing information',
            lambda: utils.create_input_files.create_sample_processing_input(
                metadata_filepath=clean_metadata_filepath,
                exclude_sampleids_list=exclude_sampleids_list,
                flowcells_directory=flowcells_directory, 
                genome_shorthand=genome_shorthand,
                dtype_mapping=dtype_mapping,
                output_filepath=sample_processing_input_filepath
                )
        )

        run_step(
            'Run to create a bed file dictating the variant positions matching the STITCH pos file',
            lambda: utils.create_samtools_depth_bed(
                positions_directory=positions_directory,
                output_directory=f"{round_directory}/output/03_sample_processing"
                )
        )

##############################################################################################################################################################################

if __name__ == "__main__":
    main()