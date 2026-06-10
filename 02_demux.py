#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 2: Demultiplex library-level FASTQs to sample-level FASTQs
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Library-level FASTQs
    - Demux input file (generated from step 1)
    - Demux array file (generated from step 1)
    - Library sheet (generated from step 1)
    - FASTQ files (Illumina paired-end short reads)

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Sample-level FASTQs
    - Library-level metrics file

Assumptions:
    - Read structures

"""

##############################################################################################################################################################################

import user
import subprocess
import pandas as pd
import sys
from time import time
from utils.print_aesthetics import *

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)

# === Load required modules (CPU, GCC, Java) for fgbio ===
subprocess.run("module load cpu/0.17.3 && module load gcc/10.2.0-2ml3m2l && module load openjdk/1.8.0_265-b01-7lpdkk2", shell=True)

# === Load system arguments ===
array_taskid = int(sys.argv[1]) 
ncpu = int(sys.argv[2])
mem_per_cpu = int(sys.argv[3])
totalmem = int(ncpu * mem_per_cpu)
javamem = f"{totalmem}G"

print(f"array_taskid : {array_taskid}")
print(f"ncpu : {ncpu}")
print(f"mem_per_cpu : {mem_per_cpu}")
print(f"totalmem : {totalmem}")
print(f"javamem : {javamem}")

# === Declare data type of special columns ===
dtype_mapping={'rfid':str}

# === Extract info for the specific library being processed ===
input_filepath = f"{round_directory}/output/02_demux/{current_round}_demux_input.csv"
demux_input = pd.read_csv(input_filepath, dtype=dtype_mapping).reset_index(drop=True)
row = demux_input.loc[demux_input['idx'] == array_taskid].iloc[0]
sequencing_method = row['seq_method']
library_sheet_filepath = row['library_sheet_filepath']
fq1_filepath = row['fq1_filepath']
fq2_filepath = row['fq2_filepath']
library_demux_outdir = row['library_demux_outdir']
library_demux_metrics_filepath = row['library_demux_metrics_filepath']

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
    print(elapsed)

def sh(cmd):
    """
    cmd: list[str]
    """
    subprocess.run(cmd, check=True)

##############################################################################################################################################################################

def riptide_block():
    run_step(
        "Demultiplex riptide library-level fastqs into sample-level fastqs",
        lambda: sh([
        "java", f"-Xmx{javamem}", "-XX:+AggressiveOpts", "-XX:+AggressiveHeap",
        "-jar", fgbio, "DemuxFastqs",
        "--inputs", fq1_filepath, fq2_filepath,
        "--metadata", library_sheet_filepath,
        "--read-structures", "8B12M+T", "8M+T",
        "--output-type=Fastq",
        "--threads", str(ncpu),
        "--output", library_demux_outdir,
        "--metrics", library_demux_metrics_filepath,
        ])
    )


def gbs_block():
    run_step(
        "Demultiplex gbs library-level fastqs into sample-level fastqs",
        lambda: print(
        """
        All ddGBS samples have been demultiplexed. The original flowcell-id level fastq files have been deleted.
        There is no need to run fastx. Below is the legacy code for fastx: 
        zcat $fastq_dir/$fastq_file | 
        $fastx_barcode_splitter --bcfile $barcode_file --prefix $demux_dir/fastq/ --bol --suffix .fq
        """
        )
    )


def main():

    if sequencing_method in ['riptide', 'LSW_lcwgs']:
        riptide_block()
    elif sequencing_method == 'gbs':
        gbs_block()
    else:
        print(f"{FAIL}: Invalid sequencing method.")
        sys.exit(1)



##############################################################################################################################################################################

if __name__ == "__main__":
    main()