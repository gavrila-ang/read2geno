#!/usr/bin/env python
import subprocess
import pandas as pd
import sys
import os
from time import time

# === Load script arguments ===
code_path = str(sys.argv[1]) 
username = str(sys.argv[2]) 
stitch_input_filepath = str(sys.argv[3]) 
stitch_arrays_filepath = str(sys.argv[4])
slurm_array_task_id = int(sys.argv[5])
slurm_job_id = int(sys.argv[6])
chunk_vcfs_output_directory = str(sys.argv[7])
temporary_directory = str(sys.argv[8]) 
nCore = int(sys.argv[9]) 
output_HD = str(sys.argv[10]) 
phased_gt = str(sys.argv[11])
bcftools_path = str(sys.argv[12])

# === Read in stitch input file ===
stitch_input_df = pd.read_csv(stitch_input_filepath).reset_index(drop=True)

# === Read in stitch arrays file ===
stitch_arrays_df = pd.read_csv(stitch_arrays_filepath).reset_index(drop=True)

# === Extract the row for the chunk being run/rerun ===
row = stitch_input_df[stitch_input_df['idx'] == slurm_array_task_id]

# === Declare STITCH arguments ===
print(f"output_filepath: {row['output_filepath'].item()}")
print(f"SLURM_ARRAY_TASK_ID: {slurm_array_task_id}")
print(f"chr: {row['chr'].item()}")
print(f"K: {row['k'].item()}")
print(f"nGen: {row['nGen'].item()}")
print(f"niterations: {row['niterations'].item()}")
print(f"method: {row['method'].item()}")
print(f"output_directory: {row['output_directory'].item()}")
print(f"temporary_directory: {temporary_directory}")
print(f"bamlist: {row['bamlist'].item()}")
print(f"sampleNames_file: {row['sampleNames_file'].item()}")
print(f"posfile: {row['posfile'].item()}")
print(f"regionStart: {row['regionStart'].item()}")
print(f"regionEnd: {row['regionEnd'].item()}")
print(f"nCore: {nCore}")
print(f"refPnls: {row['reference_panels'].item()}")
print(f"output_HD: {output_HD}")
print(f"phased_gt: {phased_gt}")

# === Ensure the output directory exists ===
os.makedirs(row['output_directory'].item(), exist_ok=True)

# === Format temporary_directory string ===
if temporary_directory == "default":
    temporary_directory = f"/scratch/{username}/job_{slurm_job_id}"


# === Find out STITCH version ===
result = subprocess.run(["conda", "list"], capture_output=True, text=True)
for line in result.stdout.splitlines():
    if "stitch" in line:
        print(line)

##############################################################################################################################################################################

# === Start Time ===
start = time()

# === Check if the output filepath exists, if not run STITCH ===
if not os.path.exists(row['output_filepath'].item()):
    print("Imputing genotypes...")
    subprocess.run([
            f"Rscript",
            f"{code_path}/utils/STITCH.R",
            f"{row['chr'].item()}",
            f"{row['k'].item()}",
            f"{row['nGen'].item()}",
            f"{row['niterations'].item()}",
            f"{row['method'].item()}",
            f"{row['output_directory'].item()}",
            f"{temporary_directory}",
            f"{row['bamlist'].item()}",
            f"{row['sampleNames_file'].item()}",
            f"{row['posfile'].item()}",
            f"{row['regionStart'].item()}",
            f"{row['regionEnd'].item()}",
            f"{nCore}",
            f"{row['reference_panels'].item()}",
            f"{output_HD}",
            f"{phased_gt}"
        ], check=True)
    print("Imputation complete.")

else:
    print("WARNING: Imputation file already exists.")

# === End Time ===
elapsed = time() - start
print(f"Time taken for imputation: {elapsed}")

##############################################################################################################################################################################

# === Start Time ===
start = time()

# === Check if the index filepath exists, if not index ===
if not os.path.exists(f"{row['output_filepath'].item()}.tbi"):
    print("Creating index file...")
    subprocess.run([
        bcftools_path, "index",
        "-t",
        "--threads", str(nCore),
        row['output_filepath'].item()
    ], check=True)
    print("Indexing complete.")

else:
    print("Indexing complete.")

# === End Time ===
elapsed = time() - start
print(f"Time taken for indexing: {elapsed}")

##############################################################################################################################################################################

# === Print out time and memory ===
print("Print out time and memory")
print("")
subprocess.run(['seff', str(slurm_job_id)], text=True)
print("")
subprocess.run([
    'sacct', 
    '-j', str(slurm_job_id), 
    '--format=JobID,JobName,MaxRSS,ReqMem,Elapsed,State'], 
    text=True)