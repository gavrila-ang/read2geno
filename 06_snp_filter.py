#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 6: Chromosome chunk combination, PLINK conversion, and SNP filtering
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Chromosome chunk STITCH file (.vcf.gz)
    - thresholds
    - dscordant snps 

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Chromosome-level VCF files (.vcf.gz)
    - Chromosome-level PLINK files 
    - List of bad SNPs
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
stitch_directory = str(sys.argv[2]) 
outdir = str(sys.argv[3])

print(f"chromvar: {chromvar}")
print(f"chunk_outdir: {stitch_directory}/stitch_chrom_chunks_vcfs")
print(f"outdir: {outdir}")

# === Declare input/output filepaths ===
chunk_outdir = f"{stitch_directory}/stitch_chrom_chunks_vcfs"
clean_metadata_filepath = f"{round_directory}/output/01_metqc/flowcell_sheets/{current_round}_metadata_cleaned.csv"

# === Declare data type of special columns ===
dtype_mapping={'rfid':str}

# === Define main filename prefixes ===
raw_prefix = f"{outdir}/chr{chromvar}_allsexes/{current_round}_chr{chromvar}_allsexes_raw"
snp_filtered_prefix = f"{outdir}/chr{chromvar}_allsexes/{current_round}_chr{chromvar}_allsexes_snp_filtered"
print(f"raw_prefix: {raw_prefix}")
print(f"snp_filtered_prefix: {snp_filtered_prefix}")

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

def index_vcf(bcf_filepath):
    # Index a bcf file
    run_step(
        f'Index {bcf_filepath}',
        lambda: sh([
            bcftools, "index", 
            "-c", 
            "--threads", ncpu, 
            bcf_filepath
        ])
    )


def bcf2ped(input_bcf=f"{raw_prefix}.bcf", plink_prefix=raw_prefix):

    # === Command to convert VCF file into PLINK (.pgen) ===
    run_step(
        f'Command to convert VCF file into PLINK (.pgen)',
        lambda: sh([
            plink2,
            "--bcf", input_bcf,
            "--thread-num", ncpu,
            # "--set-missing-var-ids", "@:#",
            "--set-all-var-ids", "@:#",            
            "--make-pgen",
            "--update-sex", f"{round_directory}/output/06_snp_filter_HD/{current_round}_plink_allsexes_sex_data.csv",
            "--keep-allele-order",
            "--output-chr", "chrM",
            "--chr-set", "20", "no-xy",
            "--out", plink_prefix,
        ])
    )

def filter_snps(raw_prefix=raw_prefix, snp_filtered_prefix=snp_filtered_prefix, bad_snps=f"{raw_prefix}_bad_snps.tsv"):

    # == Command to filter SNPs from a PLINK (.pgen) file
    run_step(
        f'Command to filter SNPs from a PLINK (.pgen) file',
        lambda: sh([
            plink2,
            "--pfile", raw_prefix,
            "--thread-num", ncpu,
            "--exclude", bad_snps,
            "--make-pgen",
            "--out", snp_filtered_prefix
        ])
    )


##############################################################################################################################################################################    

# === Blocks ===

def autosome_block():
    
    # === Combine chunks of autosomal and mitochondrial chromosomes ===
    tmpstr = glob.glob(f"{chunk_outdir}/chr{chromvar}_allsexes/stitch.chr{chromvar}.*.vcf.gz")    
    run_step(
        f'Combine chunks of chromosome {chromvar}',
        lambda: sh([
            bcftools, "concat", 
            "--threads", ncpu,
            "--no-version",
            "-a", 
            "-d", "none",
            "-Ob", 
            "-o", f"{raw_prefix}.bcf",
        ] + tmpstr)
    )

    # === Index the file ===
    index_vcf(f"{raw_prefix}.bcf")


def chromX_block():

    # === Make the output subdirectories ===
    os.makedirs(f"{outdir}/chrX_female", exist_ok = True)
    os.makedirs(f"{outdir}/chrX_male", exist_ok = True)
    os.makedirs(f"{outdir}/chrX_allsexes", exist_ok = True)

    # === Combine chunks of chrX containing only female samples ===
    tmpstr = glob.glob(f"{chunk_outdir}/chrX_female/stitch.chrX.*.vcf.gz")
    run_step(
        f'Combine chunks of chromosome {chromvar}, sex female',
        lambda: sh([
            bcftools, "concat", 
            "--threads", ncpu,
            "--no-version",
            "-a", 
            "-d", "none",
            "-Ob", 
            "-o", f"{outdir}/chrX_female/{current_round}_chrX_female_raw.bcf",
        ] + tmpstr)
    )
    # === Index the file ===
    index_vcf(f"{outdir}/chrX_female/{current_round}_chrX_female_raw.bcf")
    
    # === Combine chunks of chrX containing only male samples ===
    tmpstr = glob.glob(f"{chunk_outdir}/chrX_male/stitch.chrX.*.vcf.gz")    
    run_step(
        f'Combine chunks of chromosome {chromvar}, sex male',
        lambda: sh([
            bcftools, "concat", 
            "--threads", ncpu,
            "--no-version",
            "-a", 
            "-d", "none",
            "-Ob", 
            "-o", f"{outdir}/chrX_male/{current_round}_chrX_male_raw.bcf",
        ] + tmpstr)
    )
    # === Index the file ===
    index_vcf(f"{outdir}/chrX_male/{current_round}_chrX_male_raw.bcf")

    # === Merge chrX of females and males ===
    # === Allow for multiallelic SNP records ===
    # === Keep SNP INFO records of both males and females ===
    run_step(
        f'Combine chrX of males and females',
        lambda: sh([
            bcftools, "merge",
            "--threads", ncpu,
            "--merge", "snps",
            "--info-rules", "EAF:join,INFO_SCORE:join,HWE:join,ERC:join,EAC:join,PAF:join,REF_PANEL:join",
            f"{outdir}/chrX_female/{current_round}_chrX_female_raw.bcf",
            f"{outdir}/chrX_male/{current_round}_chrX_male_raw.bcf",
            "-Ob",
            "-o", f"{raw_prefix}_sampleid_unordered.bcf"
        ])
    )    

    # === Index the file ===
    index_vcf(f"{raw_prefix}_sampleid_unordered.bcf")
    
    # === Reorder samples in the chrX VCF file to match the rest of the chromosomes ===
    run_step(
        f'Reorder samples in chrX VCF file',
        lambda: sh([
            bcftools, "view",
            "--threads", ncpu,
            "-S", f"{stitch_directory}/stitch_lists/{current_round}_allsexes_sampleid_list.txt",
            "-Ob", 
            "-o", f"{raw_prefix}.bcf",
            f"{raw_prefix}_sampleid_unordered.bcf"
        ])
    )

    # === Index the file ===
    index_vcf(f"{raw_prefix}.bcf")


def chromY_block(exclude_female_samples=True):

    if exclude_female_samples:

        print("WARNING: Excluding female samples, will be unable to concatenate with chrY the rest of the chromosomes.")

        # === Make the output subdirectories ===
        os.makedirs(f"{outdir}/chrY_allsexes", exist_ok = True)        

        # === Combine chunks of chrY containing only male samples, allow for no duplicate SNPs ===
        tmpstr = glob.glob(f"{chunk_outdir}/chrY_male/stitch.chrY.*.vcf.gz")    
        run_step(
            f'Combine chunks of chr{chromvar}, sex male',
            lambda: sh([
                bcftools, "concat", 
                "--threads", ncpu,
                "--no-version",
                "-a", 
                "-d", "none",
                "-Ob", 
                "-o", f"{raw_prefix}.bcf",
            ] + tmpstr)
        )
        # === Index the file ===
        index_vcf(f"{raw_prefix}.bcf")

    else:
        # === Make the output subdirectories ===
        os.makedirs(f"{outdir}/chrY_male", exist_ok = True)    
        os.makedirs(f"{outdir}/chrY_female", exist_ok = True)        
        os.makedirs(f"{outdir}/chrY_allsexes", exist_ok = True)


        # === Combine chunks of chrY containing only male samples, allow for no duplicate SNPs ===
        tmpstr = glob.glob(f"{chunk_outdir}/chrY_male/stitch.chrY.*.vcf.gz")    
        run_step(
            f'Combine chunks of chr{chromvar}, sex male',
            lambda: sh([
                bcftools, "concat", 
                "--threads", ncpu,
                "--no-version",
                "-a", 
                "-d", "none",
                "-Ob", 
                "-o", f"{outdir}/chrY_male/{current_round}_chrY_male_raw.bcf",
            ] + tmpstr)
        )
        # === Index the file ===
        index_vcf(f"{outdir}/chrY_male/{current_round}_chrY_male_raw.bcf")


        # === Create snp list from chrY male unfiltered VCF' ===

        run_step(
            'Create snp list from chrY male unfiltered VCF',
            lambda: sh([
                bcftools, "query", "-f", "'%CHROM\t%POS\t%ID\t%REF\t%ALT\n'", 
                f"{outdir}/chrY_male/{current_round}_chrY_male_raw.bcf"
            ], stdout_path=f"{outdir}/chrY_female/snp_list.tsv")
        )

        # === Create an empty chrY VCF file containing only females ===
        run_step(
            f'Create an empty chrY VCF file containing only females',
            lambda: sh([
                f"{code}/utils/create_empty_vcf_with_info.sh",
                f"{outdir}/chrY_female/snp_list.tsv", # snp_list, 
                f"{round_directory}/output/05_stitch/stitch_lists/{current_round}_female_sampleid_list.txt", #sample_list, 
                f"{outdir}/chrY_female/{current_round}_chrY_female_raw.vcf", # outvcf,
                bcftools
            ])
        )

        # === Convert female chrY VCF file into a BCF file
        run_step(
            f'Convert female chrY VCF file into a BCF file',
            lambda: sh([
                bcftools, "view",
                "-Ob",
                "-o" f"{outdir}/chrY_female/{current_round}_chrY_female_raw.vcf"
            ])
        )

        # === Index the file ===
        index_vcf(f"{outdir}/chrY_female/{current_round}_chrY_female_raw.bcf")

        # === Merge chrY of females and males ===
        # === Allow for multiallelic SNP records ===
        # === --info-rules not set, so info fielfds with no specified rule will take the value ===
        # === from the first input file (chrY male vcf) ===
        run_step(
            f'Combine chrY of males and females',
            lambda: sh([
                bcftools, "merge",
                "--threads", ncpu,
                "--merge", "snps",
                f"{outdir}/chrY_male/{current_round}_chrY_male_raw.bcf",
                f"{outdir}/chrY_female/{current_round}_chrY_female_raw.bcf",
                "-Ob", 
                "-o", f"{raw_prefix}_sampleid_unordered.bcf",
            ])
        )    

        # === Index the file ===
        index_vcf(f"{raw_prefix}_sampleid_unordered.bcf")

        # === Reorder samples in the chrY VCF file to match the rest of the chromosomes ===
        run_step(
            f'Reorder samples in chrY VCF file',
            lambda: sh([
                bcftools, "view",
                "--threads", ncpu,
                "-S", f"{round_directory}/output/05_stitch/stitch_lists/{current_round}_allsexes_sampleid_list.txt",
                "-Ob", 
                "-o", f"{raw_prefix}.bcf",
                f"{raw_prefix}_sampleid_unordered.bcf",
            ])
        )

        # === Index the file ===
        index_vcf(f"{raw_prefix}.bcf")



def snpqc_block(input_bcf=f"{raw_prefix}.bcf", plink_prefix=raw_prefix, run_info=True):

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


def sampleqc_block(plink_prefix=raw_prefix, chromvar=chromvar):

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


##############################################################################################################################################################################
            

def main():

    print(f"Running SNP filtering on {chromvar}, allsexes....")    

    # === Make output subdirectories ===
    os.makedirs(f"{outdir}/chr{chromvar}_allsexes", exist_ok = True)    

    # === VCF chunk merging is chromosome-specific ===
    if chromvar not in ['X', 'Y']:
        autosome_block()
    elif chromvar == 'X':
        chromX_block()
    elif chromvar == 'Y':
        chromY_block()
    else:
        print(f"Invalid chromvar : {chromvar}")

    === Calculate quality control statistics for snps and sample from raw STITCH output ===

    # Create raw stats
    bcf2ped(input_bcf=f"{raw_prefix}.bcf", plink_prefix=raw_prefix)
    snpqc_block(input_bcf=f"{raw_prefix}.bcf", plink_prefix=raw_prefix, run_info=True)
    sampleqc_block(plink_prefix=raw_prefix, chromvar=chromvar)
    utils.snpqc_wrangling.bad_snps(plink_prefix=raw_prefix, novel_unknown_snps_filepath=novel_unknown_snps_filepath, chromvar=chromvar)

    # Create snp filtered stats
    filter_snps(raw_prefix=raw_prefix, snp_filtered_prefix=snp_filtered_prefix, bad_snps=f"{raw_prefix}_bad_snps.tsv")
    sampleqc_block(plink_prefix=snp_filtered_prefix, chromvar=chromvar)


    # === Print finishing message ===
    print(f"Completed merging the chunks of and filtering chromosome {chromvar}.")

##############################################################################################################################################################################

if __name__ == "__main__":
    main()



    
