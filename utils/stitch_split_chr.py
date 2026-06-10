#!/usr/bin/env python
import subprocess
import sys


def run_stitch_split_chr(code_path, positions_directory, output_directory):

    """
    A Python-function that runs the STITCH_split_chr.R script that splits chromosomes into chunks for STITCH arrays.
    Iterates through all the chromosomes and creates file for each. 
    
    Parameters:
    - code_path (str): Path to the code stated in the `pipeline_args.sh` file.
    - positions_directory (str): Path to the STITCH positions file stated in the `pipeline_args.sh` file 
    - output_directory (str): Path to where the output files will be saved.
    """
    
    # List of chromosomes in the reference genome
    chrom_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                    'X', 'Y', 'M'
                ]

    # Create positions where the chromosome will be divide to create chunks
    # The starting position is where the first SNP in chromosome positions file is found
    # The ending position is where the last SNP in chromosome positions file is found
    for chrom in chrom_list:
        subprocess.run(
            f"Rscript {code_path}/STITCH_split_chr.R {positions_directory}/chr{chrom}_STITCH_baud_pos {output_directory}/chr{chrom}_chunks.txt",
            shell=True, check=True
        )
        
        print(f"Finished splitting chromosome : {chrom}!")
