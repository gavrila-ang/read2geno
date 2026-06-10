#!/usr/bin/env python3

""" 
Palmer Lab Genotyping and Validation Pipeline
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Metadata file (obtained from https://irs.ratgenes.org)
    - FASTQ files (Illumina paired-end short reads)
    - Reference genome (.fa or .fasta)
    - SNP Positions File (tab-delimited file)
    - Imputation Reference Panel (.hap, .legend, .samples)

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Variant Call Format files (.vcf)
    - PLINK files (.bim, .fam, .bed)

Assumptions:
    - Consistent sequencing protocol across samples
    - Custom input file directory structure

"""

##############################################################################################################################################################################

# === Load python packages ===
import argparse
import os
import subprocess
import shlex
from datetime import datetime
import user
import utils
import sys
import math
import pandas as pd
from tqdm import tqdm

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)

# === Load current date and time ===
dtnow = datetime.now().strftime("%Y%m%d_%H%M%S")

##############################################################################################################################################################################

# === Function to parse command-line arguments to run specific pipeline steps ===
def parse_args(args):

    parser = argparse.ArgumentParser(description="Select which pipeline step to run.")

    parser.add_argument(
        "--step",
        choices=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
        required=True,
        help="python pipeline.py --step 1"
    )

    return parser.parse_args(args)

##############################################################################################################################################################################

# === Load pipeline configurations ===
def main():

    # === Print welcome message ===
    print(f"Welcome to the Palmer Lab Genotyping and Validation Pipeline!")
    print(f"Est. 2021")
    print("")

    # === Load command-line arguments ===
    args = parse_args(sys.argv[1:]) #[0] is pipeline.py


    # === Run Step 1: Quality control of metadata file ===
    if args.step == "1":
        subprocess.run(["conda", "run", "-n", conda_env, "python", f"{code}/01_metqc.py"], check=True)

    # === Run Step 2: Demultiplex library-level fastqs to sample-level fastqs ===
    elif args.step == "2":

        # === Check if demultiplexing arrays finished successfully ===
        subprocess.run(["conda", "run", "-n", conda_env, "python", f"{code}/02_demux_check.py"], check=True)

        # === Create a python list from the array file created ===
        array_filepath = f"{round_directory}/output/02_demux/02_demux_arrays.csv"
        tempdf = pd.read_csv(array_filepath, header=None, names=['array_val'], dtype={'array_val': str})
        array_list = tempdf['array_val'].tolist()
        array_str = ','.join(array_list)

        # === Execute parallelised library demultiplexing tasks ===
        if len(array_list) != 0:       
            print(f"Running Step 2.2: Library demultiplexing for arrays {array_str}....")              
            for array_val in array_list:
                subprocess.run([
                    "sbatch", "-J", f"02_demux_arr{array_val}",
                    "-o", f"{round_directory}/logs/02_demux/{dtnow}-%A-%a.o",
                    "-e", f"{round_directory}/logs/02_demux/{dtnow}-%A-%a.e",
                    "-N", "1", "-c", str(ncpu), f"--mem-per-cpu={mpc}G",
                    "-t", "24:00:00",
                    "-A", allocation, "-p", partition, "-q", qos,
                    f"--exclude={exclude_nodes}",
                    "--wrap", shlex.join(["python", f"{code}/02_demux.py", str(array_val), str(ncpu), str(mpc)])
                    ], check=True)

    # === Run Step 3: Clean sample-level fastqs, align to reference genome and test alignment quality ===
    elif args.step == "3":

        # === Check if aligning arrays finished successfully ===
        subprocess.run(["conda", "run", "-n", conda_env, "python", f"{code}/03_sample_processing_check.py"], check=True)

        # === Create a python list from the array file created ===
        array_filepath = f"{round_directory}/output/03_sample_processing/03_sample_processing_arrays.csv"
        tempdf = pd.read_csv(array_filepath, dtype={'idx': str})
        array_list = tempdf.head(998)['idx'].tolist() # tscc maximum 999 arrays
        array_str = ','.join(array_list)

        # === Execute parallelised sample processing tasks ===
        if len(array_list) > 0:
            print(f"Running Step 3: Sample processing for {len(array_list)} samples")
            for idx in array_list:
                subprocess.run([
                    "sbatch", "-J", f"03_sample_processing_arr{idx}",
                    "-o", f"{round_directory}/logs/03_sample_processing/{dtnow}-%A-%a.o",
                    "-e", f"{round_directory}/logs/03_sample_processing/{dtnow}-%A-%a.e",
                    "-N", "1", "-c", str(ncpu), f"--mem-per-cpu={mpc}G",
                    "-t", "12:00:00",
                    "-A", allocation, "-p", partition, "-q", qos,
                    f"--exclude={exclude_nodes}",
                    "--wrap", f"bash -c 'source ~/.bashrc && conda activate {conda_env} && python {code}/03_sample_processing.py {idx} {ncpu} {mpc}'"
                    ], check=True)

    elif args.step == "4":
        run_type = "abv"
        gbs_sexqc_limit = 0.95
        riptide_sexqc_limit = 0.95
        skip_sexqc = True
        stitch_directory = f"{round_directory}/output/05_stitch"
        print(f"Running Step 4: Alignment QC....")
        subprocess.run([
                "sbatch", "-J", "04_alignqc",
                "-o", f"{round_directory}/logs/04_align_qc/{dtnow}-%A-%a.o",
                "-e", f"{round_directory}/logs/04_align_qc/{dtnow}-%A-%a.e",
                "-N", "1", "-c", str(ncpu), f"--mem-per-cpu={mpc}G",
                "-t", "1:00:00",
                "-A", allocation, "-p", partition, "-q", qos,
                f"--exclude={exclude_nodes}",
                "--wrap", f"bash -c 'source ~/.bashrc && conda activate {conda_env} && python {code}/04_alignqc.py {run_type} {gbs_sexqc_limit} {riptide_sexqc_limit} {skip_sexqc} {stitch_directory}'"
                ], check=True)

    elif args.step == "5":

        # STITCH user arguments
        chunk_vcfs_output_directory = f"{round_directory}/output/05_stitch/stitch_chrom_chunks_vcfs"
        temporary_directory = "default"  
        output_hd = "False"
        phased_gt = "False"
        stitch_directory = f"{round_directory}/output/05_stitch"

        print(f"Running Step 5: STITCH....")
        subprocess.run(
            f"python {code}/05_stitch.py "
            f"{chunk_vcfs_output_directory} "
            f"{temporary_directory} "
            f"{output_hd} "
            f"{stitch_directory} "
            f"{phased_gt}",
            shell=True
         )


    elif args.step == "6":

        # Additional user arguments
        stitch_directory = f"{round_directory}/output/05_stitch"
        outdir = f"{round_directory}/output/06_snp_filter"
        temporary_directory = f"{round_directory}/output/tempdir"
        pedigree_filepath = f"{round_directory}/output/06_snp_filter/{current_round}_pedigree.csv"

        chrom_list = list(range(1, 21)) + ['X', 'Y', 'M']                
        # chrom_list = [12, 'X', 'Y', 'M']

        for chromvar in chrom_list:
            print(f"Running Step 6: SNP Filtering on {chromvar}....")
            subprocess.run([
                "sbatch", "-J", f"06_snpqc_chr{chromvar}",
                "-o", f"{round_directory}/logs/06_snp_filter/{dtnow}-%A-%a.o",
                "-e", f"{round_directory}/logs/06_snp_filter/{dtnow}-%A-%a.e",
                "-N", "1", "-c", ncpu, f"--mem-per-cpu={mpc}G", 
                "-t", "12:00:00",
                "-A", allocation, "-p", partition, "-q", qos,
                f"--exclude={exclude_nodes}",
                "--wrap", f"bash -c 'source ~/.bashrc && conda activate {conda_env} && python {code}/06_snp_filter.py {chromvar} {stitch_directory} {outdir}'"
                ], check=True)


    elif args.step == "7":
        
        # Additional user arguments
        snp_filtered_outdir = f"{round_directory}/output/06_snp_filter"
        snp_and_sample_filtered_outdir = f"{round_directory}/output/07_sample_filter"
        mer_threshold = 0.05    
        missingness_threshold = 10    

        # Short python function to combine the create_sampleqc_summary across all chromosomes to create an aggregated bad_samples.tsv
        utils.sampleqc_wrangling.bad_samples(snp_filtered_outdir=snp_filtered_outdir, snp_and_sample_filtered_outdir=snp_and_sample_filtered_outdir, current_round=current_round)

        chrom_list = list(range(1, 21)) + ['X', 'Y', 'M']                

        for chromvar in chrom_list:
            print(f"Running Step 7: Sample Filtering on {chromvar}....")
            subprocess.run([
                "sbatch", "-J", f"07_sampleqc_chr{chromvar}",
                "-o", f"{round_directory}/logs/07_sample_filter/{dtnow}-%A-%a.o",
                "-e", f"{round_directory}/logs/07_sample_filter/{dtnow}-%A-%a.e",
                "-N", "1", "-c", ncpu, f"--mem-per-cpu={mpc}G", 
                "-t", "12:00:00",
                "-A", allocation, "-p", partition, "-q", qos,
                f"--exclude={exclude_nodes}",
                "--wrap", f"bash -c 'source ~/.bashrc && conda activate {conda_env} && python {code}/07_sample_filter.py {chromvar} {snp_filtered_outdir} {snp_and_sample_filtered_outdir}'"
                ], check=True)


    # === Run Step 8: Generate database-ready QC log files ===
    elif args.step == "8":

        # Additional user arguments
        schema_output_directory = f"{round_directory}/output/08_schema"
        genome_size = 2_845_466_930  # GRCr8 rat reference genome size (bp)

        os.makedirs(f"{round_directory}/logs/08_schema", exist_ok=True)
        os.makedirs(schema_output_directory, exist_ok=True)

        print(f"Running Step 8: Generating database QC logs....")
        subprocess.run([
            "sbatch", "-J", "08_schema",
            "-o", f"{round_directory}/logs/08_schema/{dtnow}-%A-%a.o",
            "-e", f"{round_directory}/logs/08_schema/{dtnow}-%A-%a.e",
            "-N", "1", "-c", str(ncpu), f"--mem-per-cpu={mpc}G",
            "-t", "12:00:00",
            "-A", allocation, "-p", partition, "-q", qos,
            f"--exclude={exclude_nodes}",
            "--wrap", f"bash -c 'source ~/.bashrc && conda activate {conda_env} && python {code}/08_schema.py {schema_output_directory} {genome_size}'"
            ], check=True)


    # === Run Step 9: Generate publication-ready QC figures ===
    elif args.step == "9":

        # Additional user arguments
        database_logs_directory = f"{round_directory}/output/08_schema"
        report_output_directory = f"{round_directory}/output/09_publication"

        os.makedirs(f"{round_directory}/logs/09_publication", exist_ok=True)
        os.makedirs(report_output_directory, exist_ok=True)

        print(f"Running Step 9: Generating publication figures....")
        subprocess.run([
            "sbatch", "-J", "09_publication",
            "-o", f"{round_directory}/logs/09_publication/{dtnow}-%A-%a.o",
            "-e", f"{round_directory}/logs/09_publication/{dtnow}-%A-%a.e",
            "-N", "1", "-c", str(ncpu), f"--mem-per-cpu={mpc}G",
            "-t", "12:00:00",
            "-A", allocation, "-p", partition, "-q", qos,
            f"--exclude={exclude_nodes}",
            "--wrap", f"bash -c 'source ~/.bashrc && conda activate {conda_env} && python {code}/09_publication.py {database_logs_directory} {report_output_directory}'"
            ], check=True)



##############################################################################################################################################################################


if __name__ == '__main__':
    main()
