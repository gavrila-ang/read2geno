#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 5: Genotyping by imputation using STITCH
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Alignemnt files (.bam)
    - Bam list (.csv)
    - Sample list (.csv)
    - STITCH input file
    - STITCH array file
    - Reference genome

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - VCF files (.vcf.gz)

Assumptions:
    - Custom input file directory structure

"""

##############################################################################################################################################################################

# === Load python packages ===
import user
import utils
import subprocess
import pandas as pd
import sys
from datetime import datetime


# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Load system arguments ===
chunk_vcfs_output_directory = str(sys.argv[1]) 
temporary_directory = str(sys.argv[2]) 
output_HD = str(sys.argv[3])
stitch_directory = str(sys.argv[4])
phased_gt = str(sys.argv[5])

print(f"chunk_vcfs_output_directory: {chunk_vcfs_output_directory}")
print(f"temporary_directory: {temporary_directory}")
print(f"output_HD: {output_HD}")
print(f"stitch_directory: {stitch_directory}")
print(f"phased_gt: {phased_gt}")

# === Declare input/output filepaths ===
clean_metadata_filepath = f"{round_directory}/output/01_metqc/flowcell_sheets/{current_round}_metadata_cleaned.csv"
stitch_input_filepath = f"{stitch_directory}/{current_round}_stitch_input.csv"
stitch_arrays_filepath = f"{stitch_directory}/05_stitch_arrays.csv"

# === Declare data type of special columns ===
dtype_mapping={'rfid':str}

##############################################################################################################################################################################

# === For each row in stitch input, has the output file been produced; produce rerun array file ===
utils.create_stitch_arrays(metadata_filepath=clean_metadata_filepath, 
                            round_directory=round_directory,
                            current_round=current_round,
                            stitch_input_filepath=stitch_input_filepath, 
                            output_filepath=stitch_arrays_filepath, 
                            exclude_sampleids_list=exclude_sampleids_list, 
                            dtype_mapping=dtype_mapping)

# === Loop through user input range of rerun array file ===
stitch_arrays = pd.read_csv(stitch_arrays_filepath, dtype={'idx':'Int64', 'chr':'str', 'sex':'str'}).reset_index()
array_list = stitch_arrays['idx'].values.tolist()

# === Genotype the chromosome chunks! ====
if array_list:
    print(f"\nSubmitting arrays for STITCH....")
for array in array_list:
    dtnow = datetime.now().strftime("%Y%m%d_%H%M%S")
    subprocess.run([
            "sbatch",
            "-J", f"05_stitch_{array}",
            "-o", f"{round_directory}/logs/05_stitch/{dtnow}-%A-%a.o",
            "-e", f"{round_directory}/logs/05_stitch/{dtnow}-%A-%a.e",
            "-N", "1",
            "--array", str(array),
            "-c", str(ncpu),
            f"--mem-per-cpu={mpc}G",
            "-t", "8:00:00",
            "-A", allocation,
            "-p", partition,
            "-q", qos,
            "--exclude", exclude_nodes,
            "--wrap", (
                f"source $HOME/miniconda3/etc/profile.d/conda.sh && "
                f"conda activate {stitch_env} && "
                f"python {code}/utils/STITCH.py "
                f"{code} {username} {stitch_input_filepath} {stitch_arrays_filepath} "
                f"$((SLURM_ARRAY_TASK_ID)) $((SLURM_JOBID)) {chunk_vcfs_output_directory} {temporary_directory} "
                f"{ncpu} {output_HD} {phased_gt} {bcftools}"
            )
        ],
        stdout=subprocess.PIPE,
        text=True
    )
