#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 1: Quality control of metadata file
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Metadata file (obtained from https://irs.ratgenes.org)
    - FASTQ files (Illumina paired-end short reads)

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Quality control files
    - Cleaned metadata file

Assumptions:
    - Custom input file directory structure

"""

##############################################################################################################################################################################

import user
import utils
import subprocess
from time import time
from utils.print_aesthetics import *

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Declare output directories ===
html_fig_dir = f"{round_directory}/output/01_metqc/html_figures"
log_dir = f"{round_directory}/logs/01_metqc"

# === Declare main output filepath ===
clean_metadata_filepath = f"{round_directory}/output/01_metqc/flowcell_sheets/{current_round}_metadata_cleaned.csv"

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

    print(f"Running Step 1: Quality control of metadata file....")

    run_step(
        'Create log directory of current round',
        lambda: utils.create_dirs.create_round_dirs(
            round_dir=round_directory,
            code=code,
            structure_type='logs'
        )
    )

    run_step(
        'Create output directory of current round',
        lambda: utils.create_dirs.create_round_dirs(
            round_dir=round_directory,
            code=code,
            structure_type='output'
        )
    )

    run_step(
        'Create output and log directory of new flowcells in current round',
        lambda: utils.create_dirs.create_flowcell_dirs(
            metadata_filepath=initial_metadata_filepath,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            code=code,
            log_dir=log_dir,
            current_round=current_round
        )
    )

    run_step(
        'Ensure essential columns exist in the metadata in metadata',
        lambda: utils.metadata_qc.essential_columns(
            metadata_filepath=initial_metadata_filepath,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Rename columns in metadata to become fgbio-compatible',
        lambda: utils.metadata_qc.rename_columns(
            metadata_filepath=initial_metadata_filepath,
            dtype_mapping=dtype_mapping,
            output_filepath=clean_metadata_filepath
        )
    )

    run_step(
        'Create THE key column, Sample_ID',
        lambda: utils.metadata_qc.create_sampleid(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping,
            output_filepath=clean_metadata_filepath
        )
    )

    run_step(
        'Remove Sample IDs in exclude sampleids list',
        lambda: utils.metadata_qc.remove_bad_sampleids(
            metadata_filepath=clean_metadata_filepath,
            current_round=current_round,
            drop_sample_log_filepath=f"{code}/user/drop_sample_log.csv", 
            drop_library_log_filepath=f"{code}/user/drop_library_log.csv", 
            drop_project_log_filepath=f"{code}/user/drop_project_log.csv",
            reinclude_sample_log_filepath=f"{code}/user/reinclude_log.csv",
            dtype_mapping=dtype_mapping,
            output_filepath=clean_metadata_filepath
        )
    )

    run_step(
        'Find RFID with multiple trailing zeros (933000000000000)',
        lambda: utils.metadata_qc.rfid_multizero(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Find RFID in scientific notation (9.33E+14)',
        lambda: utils.metadata_qc.rfid_scinotation(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Find short RFID values in metadata',
        lambda: utils.metadata_qc.rfid_short(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping,
            log_dir=log_dir,
            current_round=current_round,
            len_limit=6
        )
    )

    run_step(
        'Check for missing RFID data in required fields',
        lambda: utils.metadata_qc.rfid_missing_data(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping,
            current_round=current_round,
            log_dir=log_dir
        )
    )

    run_step(
        'Validate FASTQ file paths',
        lambda: utils.metadata_qc.fastq_validity(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping,
            sequencing_directory=sequencing_data,
            current_round=current_round,
            log_dir=log_dir
        )
    )

    run_step(
        'Check unique library associations',
        lambda: utils.metadata_qc.library_association(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Check unique FASTQ file associations',
        lambda: utils.metadata_qc.fastq_file_association(
            metadata_filepath=clean_metadata_filepath,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Detect duplicate Sample_IDs within libraries',
        lambda: utils.metadata_qc.detect_dup_sampleids(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Create library-level flowcell sheets',
        lambda: utils.metadata_qc.make_library_sheets(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            dtype_mapping=dtype_mapping,
            output_dir=f"{round_directory}/output/01_metqc/library_sheets"
        )
    )

    run_step(
        'Create demographics summary',
        lambda: utils.metadata_plot.create_demographic_table(
            metadata_filepath=clean_metadata_filepath,
            file_prefix=f"{html_fig_dir}/{current_round}"
        )
    )

    run_step(
        'Plot demographics summary',
        lambda: utils.metadata_plot.create_demographic_image(
            current_round=current_round,
            file_prefix=f"{html_fig_dir}/{current_round}"
        )
    )

    run_step(
        'Create demultiplexing input file',
        lambda: utils.create_input_files.create_demux_input(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            dtype_mapping=dtype_mapping,
            round_directory=round_directory,
            sequencing_directory=sequencing_data,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            current_round=current_round,
            log_dir=log_dir
        )
    )

    run_step(
        'Create a file containing arrays of libraries to be demuxed',
        lambda: utils.create_rerun_arrays.create_demux_arrays(
            demux_input_filepath=f"{round_directory}/output/02_demux/{current_round}_demux_input.csv",
            metadata_filepath=clean_metadata_filepath,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            exclude_sampleids_list=exclude_sampleids_list,
            summary_filepath=f"{round_directory}/output/02_demux/undemuxed_samples.csv",
            output_filepath=f"{round_directory}/output/02_demux/02_demux_arrays.csv",
            dtype_mapping=dtype_mapping
        )
    )


    print("\nPipeline complete.")

##############################################################################################################################################################################

if __name__ == "__main__":
    main()