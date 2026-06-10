#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 7: SNP and sample filtering
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Chromosome-level PLINK files (filtered for bad SNPs)
    - List of bad SNPs
    - List of bad samples        

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Chromosome-level PLINK files 
    - SNP statistics
    - Sample statistics

Assumptions:
    - Custom input file directory structure

"""

##############################################################################################################################################################################
import user
import utils
import subprocess
import sys
import os
import pandas as pd
import glob
from datetime import datetime
from time import time


# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Load system arguments ===
chromvar = str(sys.argv[1]) 
snp_filtered_outdir = str(sys.argv[2]) 
snp_and_sample_filtered_outdir = str(sys.argv[3])

print(f"chromvar: {chromvar}")
print(f"snp_filtered_outdir: {snp_filtered_outdir}")
print(f"snp_and_sample_filtered_outdir: {snp_and_sample_filtered_outdir}")

# === Declare input/output filepaths ===
clean_metadata_filepath = f"{round_directory}/output/01_metqc/flowcell_sheets/{current_round}_metadata_cleaned.csv"

# === Declare data type of special columns ===
dtype_mapping={'rfid':str}

# === Define main filename prefixes ===
raw_prefix = f"{snp_filtered_outdir}/chr{chromvar}_allsexes/{current_round}_chr{chromvar}_allsexes_raw"
snp_filtered_prefix = f"{snp_filtered_outdir}/chr{chromvar}_allsexes/{current_round}_chr{chromvar}_allsexes_snp_filtered"
snp_and_sample_filtered_prefix = f"{snp_and_sample_filtered_outdir}/chr{chromvar}_allsexes/{current_round}_chr{chromvar}_allsexes_snp_and_sample_filtered"
print(f"snp_filtered_prefix: {snp_filtered_prefix}")
print(f"snp_and_sample_filtered_prefix: {snp_and_sample_filtered_prefix}")

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
    print(f"Time taken: {elapsed}")

def sh(cmd, stdout_path=None):
    """
    cmd: list[str]
    """    
    print(cmd)
    if stdout_path:
        with open(stdout_path, "w") as out:
            subprocess.run(cmd, check=True, stdout=out)
    else:
        subprocess.run(cmd, check=True)

##############################################################################################################################################################################

# === Functions ===

def filter_snps_and_samples(snp_filtered_prefix, snp_and_sample_filtered_prefix, bad_snps, bad_samples):

    # == Command to filter SNPs and samples from a PLINK (.pgen) file
    run_step(
        f'Command to filter SNPs and samples from a PLINK (.pgen) file',
        lambda: sh([
            plink2,
            "--pfile", snp_filtered_prefix,
            "--thread-num", ncpu,
            "--exclude", bad_snps,
            "--remove", bad_samples,
            "--make-pgen",
            "--out", snp_and_sample_filtered_prefix
        ])
    )

##############################################################################################################################################################################    

# === Blocks ===

def snpqc_block(input_bcf, plink_prefix, run_info=True):

    # === Calculate SNP INFO scores ===
    if run_info:
        run_step(
            'Calculate SNP INFO scores',
            lambda: sh([
                bcftools, "query", 
                "-HH", 
                "-f", "%CHROM\t%POS\t%CHROM:%POS\t%INFO_SCORE\n", 
                input_bcf, 
            ], stdout_path=f"{plink_prefix}.info")
        )

    # === Calculate SNP missingness rate ===
    run_step(
        'Calculate SNP missingness',
        lambda: sh([
            plink2, 
            "--thread-num", ncpu,
            "--pfile", plink_prefix,
            "--missing", "variant-only",
            "--out", plink_prefix
        ])
    )

    # === Calculate SNP minor allele frequency ===
    run_step(
        'Calculate SNP MAF ',
        lambda: sh([
            plink2, 
            "--thread-num", ncpu,
            "--pfile", plink_prefix,
            "--freq",
            "--out", plink_prefix
        ])
    )

    # === Calculate SNP hardy-weinberg equillibrium probability ===
    run_step(
        'Calculate SNP HWE',
        lambda: sh([
            plink2, 
            "--thread-num", ncpu,
            "--pfile", plink_prefix,
            "--hardy",
            "--out", plink_prefix
        ])
    )


    # # === Calculate SNP mendlian error rate ===
    # run_step(
    #     'Calculate SNP MER',
    #     lambda: sh([
    #         plink2, 
    #         "--thread-num", ncpu,
    #         "--pfile", plink_prefix,
    #         "--mendel",
    #         "--out", plink_prefix
    #     ])
    # )


def sampleqc_block(plink_prefix, chromvar):

    # === Calculate Sample heterozygosity ===
    # --het requires at least one polymorphic autosomal variant.
    if chromvar not in ['X', 'Y', 'M']:
        run_step(
            'Calculate Sample Heterozygosity',
            lambda: sh([
                plink2, 
                "--thread-num", ncpu,
                "--pfile", plink_prefix,
                "--het",
                "--out", plink_prefix
            ])
        )

   # === Calculate Sample missingness rate ===
    run_step(
        'Calculate Sample missingness',
        lambda: sh([
            plink2, 
            "--thread-num", ncpu,
            "--pfile", plink_prefix,
            "--missing", "sample-only",
            "--out", plink_prefix
        ])
    )

    # # === Calculate Sample mendlian error rate ===
    # run_step(
    #     'Calculate Sample MER',
    #     lambda: sh([
    #         plink2, 
    #         "--thread-num", ncpu,
    #         "--pfile", plink_prefix,
    #         "--mendel",
    #         "--out", plink_prefix
    #     ])
    # )


def sample_pca_block(plink_prefix):

    # === Calculate Sample PCA ===
    run_step(
        'Calculate Sample PCA',
        lambda: sh([
            plink2,
            "--thread-num", ncpu,
            "--pfile", plink_prefix,
            "--pca",
            "--out", plink_prefix
        ])
    )

    # # === Calculate Sample PCA outliers based on their sequencing methods ===
    # run_step(
    #     'Calculate sample PCA outliers',
    #     lambda: utils.sampleqc_wrangling.calculate_sample_pca_outliers(
    #         prefix=plink_prefix,
    #         metadata_filepath=clean_metadata_filepath,
    #         dtype_mapping=dtype_mapping,
    #         output_filepath=f"{plink_prefix}_outlier_samples.txt"
    #     )
    # )
    # # === Re-calculate Sample PCA without outliers ===
    # run_step(
    #     'Re-calculate Sample PCA without outliers',
    #     lambda: sh([
    #         plink2,
    #         "--thread-num", ncpu,
    #         "--pfile", plink_prefix,
    #         "--remove", f"{plink_prefix}_outlier_samples.txt",
    #         "--pca",
    #         "--out", f"{plink_prefix}_no_outlier_samples"
    #     ])
    # )

##############################################################################################################################################################################
            

def main():

    print(f"Running SNP filtering on {chromvar}, allsexes....")    

    # === Make output subdirectories ===
    os.makedirs(f"{snp_and_sample_filtered_outdir}/chr{chromvar}_allsexes", exist_ok = True)    

    # === Calculate quality control statistics for snps and sample from raw STITCH output ===

    filter_snps_and_samples(
        snp_filtered_prefix=snp_filtered_prefix, 
        snp_and_sample_filtered_prefix=snp_and_sample_filtered_prefix, 
        bad_snps=f"{raw_prefix}_bad_snps.tsv", 
        bad_samples=f"{snp_filtered_prefix}_bad_samples.tsv"
        )
    snpqc_block(input_bcf=None, plink_prefix=snp_and_sample_filtered_prefix, run_info=False)
    sampleqc_block(plink_prefix=snp_and_sample_filtered_prefix, chromvar=chromvar)
    sample_pca_block(plink_prefix=snp_and_sample_filtered_prefix)


    # === Print finishing message ===
    print(f"Completed merging the chunks of and filtering chromosome {chromvar}.")

##############################################################################################################################################################################

if __name__ == "__main__":
    main()
