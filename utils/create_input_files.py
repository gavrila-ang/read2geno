#!/usr/bin/env python
import pandas as pd
import os
import csv
import warnings
import shutil
import sys
from tqdm import tqdm
warnings.simplefilter(action='ignore', category=FutureWarning)
import subprocess
from utils.validate_paths import validate_paths
from utils.print_aesthetics import *


##############################################################################################################################################################################

def create_demux_input(metadata_filepath, exclude_sampleids_list, dtype_mapping, round_directory, sequencing_directory, flowcells_directory, genome_shorthand, current_round, log_dir):
    
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # Filter out bad samples
    metadata_df = metadata_df[~metadata_df.Sample_ID.isin(exclude_sampleids_list)]

    # Track which library_flowcell combinations have bad paths
    bad_library_flowcell = False
    
    # store demux input rows
    rows = []

    # create an empty file to store warnings
    with open(f"{log_dir}/{current_round}_demux_input_omitted_libraries.txt", 'w') as f:
        pass 

    # each demux input row pertains to a unique library+flowcell combination
    for library_flowcell, group in  metadata_df.groupby(["Library_ID", "flowcell_id"]):

        # WARNING: Assumes each Library_ID:flowcell_id combination is uniquely associated with one 'Fastq_Files' value
        # WARNING: Assumes each Library_ID:flowcell_id combintation is uniquely associated with one 'seq_method' value
        sequencing_method = group.seq_method.unique()[0]
        library_id = library_flowcell[0]
        flowcell_id = library_flowcell[1]

        if sequencing_method in ['riptide', 'LSW_lcwgs']:
            fq1_filename, fq2_filename = group.Fastq_Files.unique()[0].split(";")
        elif sequencing_method == "gbs":
            fq1_filename = group.Fastq_Files.unique().item()
            fq2_filename = "NA"
        else:
            print(f"{FAIL} ERROR: Invalid sequencing method {sequencing_method} for {library_flowcell, group}")            

        row_dict = {
            "Library_ID": library_id,
            "flowcell_id": flowcell_id,
            "seq_method": sequencing_method,
            "library_sheet_filepath": f"{round_directory}/output/01_flowcell_setup/library_sheets/{flowcell_id}_{library_id}_library_sheet.csv",
            "fq1_filepath": f"{sequencing_directory}/{flowcell_id}/{fq1_filename.strip()}",
            "fq2_filepath": f"{sequencing_directory}/{flowcell_id}/{fq2_filename.strip()}",
            "library_demux_outdir": f"{flowcells_directory}/{flowcell_id}/{genome_shorthand}/02_demux/demux/{library_id}",
            "library_demux_metrics_filepath": f"{flowcells_directory}/{flowcell_id}/{genome_shorthand}/02_demux/demux_metrics/{library_id}_demux_barcode_metrics.txt",
            "library_sheet_outdir": f"{flowcells_directory}/{flowcell_id}/{genome_shorthand}/01_flowcell_setup/library_sheets"
        }
        
        # Create library_demux_outdir
        os.makedirs(row_dict['library_demux_outdir'], exist_ok=True)
        
        # Copy library sheet from round directory to library-specific directory
        try:
            shutil.copy(row_dict['library_sheet_filepath'], row_dict['library_sheet_outdir'])
        except Exception as e:
            print(f"{FAIL} ERROR: Could not copy {row_dict['library_sheet_filepath']} to {row_dict['library_sheet_outdir']}")
            print(e)
            pass

        # Check filepaths for library
        library_sheet_output_file = f"{row_dict['library_sheet_outdir']}/{flowcell_id}_{library_id}_library_sheet.csv"
        library_demux_metrics_outdir = os.path.dirname(row_dict['library_demux_metrics_filepath'])
        filepaths_list = [row_dict['library_sheet_filepath'], 
                          row_dict['fq1_filepath'], 
                          row_dict['fq2_filepath'],
                          row_dict['library_demux_outdir'], 
                          library_demux_metrics_outdir,
                          library_sheet_output_file
                          ]

        if sequencing_method in ['riptide', 'LSW_lcwgs']:
            bad_paths = [path for path in filepaths_list if not os.path.exists(path)]
        # ddgbs samples aren't paired end sequences, so don't check fq2 path
        elif sequencing_method == "gbs":
            # don't check fq2_filepath for ddgbs samples
            bad_paths = [path for path in filepaths_list if path != 'fq2_filepath' and not os.path.exists(path)]
        else:
            print(f"{FAIL} ERROR: Invalid sequencing method {sequencing_method} for {library_flowcell, group}")

        if bad_paths:
            rows.append(row_dict)   
            bad_library_flowcell = True

            fail_message = f"{FAIL} QC FAIL: For {library_flowcell} the following paths required for demultiplexing don't exist {bad_paths}. Consider deleting from demux input file."        
            # print(fail_message)
            with open(f"{log_dir}/{current_round}_demux_input_omitted_libraries.txt", "a") as f:
                f.write(fail_message + "\n")            
            
        else:
            rows.append(row_dict)

    
    # create demux input dataframe
    demux_input = pd.DataFrame(rows)
    
    # state if all library_flowcell combinations have valid demux filepaths 
    if not bad_library_flowcell:
        print(f"{PASS} QC PASS: All library_flowcell combinations have valid demux input and output filepaths.")
    else:
        print(f"{FILE} Saving omitted libraries to {log_dir}/{current_round}_demux_input_omitted_libraries.txt")

    # save demux input file
    output_filepath = f"{round_directory}/output/02_demux/{current_round}_demux_input.csv"
    demux_input.to_csv(output_filepath, index=True, index_label="idx")
    print(f"{FILE} Saving demux input file to {output_filepath}.")



##############################################################################################################################################################################


def create_sample_processing_input(metadata_filepath, exclude_sampleids_list, dtype_mapping,
                                   flowcells_directory, genome_shorthand, output_filepath):
    """
    Creates a CSV file for downstream sample processing (cleaning, mapping), containing metadata and filepaths
    for each demultiplexed sample that should be processed.

    Parameters:
    - metadata_filepath (str): Path to the flowcell-level samplesheet CSV.
    - exclude_sampleids_list (list): Sample IDs to exclude from processing.
    - dtype_mapping (dict): Data types to use when reading the CSV.
    - flowcells_directory (str): Base directory where flowcell data is stored.
    - genome_shorthand (str): Name representing the genome reference.
    - output_filepath (str): Path to write the resulting sample processing input CSV.
    """

    # Read sample metadata
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Track which sample_ids have bad paths but aren't in exclude_sampleids_list
    bad_sample_ids = False
    
    # List to store info about each sample to process
    samples_variables = []

    for _, row in tqdm(metadata_df.iterrows(), total=len(metadata_df)):

         # List to store missing paths for a sample
        missing_paths = []

        if row.Sample_ID not in exclude_sampleids_list:

            # Construct paths for FASTQ files and output directory
            base_path = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}"
            output_dir = f"{base_path}/03_sample_processing"

            if row['seq_method'] in ['riptide', 'LSW_lcwgs']:
                fq1_filepath = f"{base_path}/02_demux/demux/{row.Library_ID}/{row.Sample_ID}-{row.Sample_Name}-{row.Sample_Barcode}_R1.fastq.gz"
                fq2_filepath = f"{base_path}/02_demux/demux/{row.Library_ID}/{row.Sample_ID}-{row.Sample_Name}-{row.Sample_Barcode}_R2.fastq.gz"                
                paths_list = [base_path, fq1_filepath, fq2_filepath, output_dir]
                missing_paths = [path for path in paths_list if not os.path.exists(path)]

            # ddgbs samples aren't paired end sequences, so don't check fq2 path
            elif row['seq_method'] == "gbs":
                fq1_filepath = f"{base_path}/02_demux/demux/{row.Library_ID}/{row.Sample_ID}.fq.gz"
                paths_list = [base_path, fq1_filepath, output_dir]
                missing_paths = [path for path in paths_list if not os.path.exists(path)]

            else:
                print(f"{FAIL} ERROR: Invalid sequencing method {row.seq_method} for {row.Sample_ID}")


            if missing_paths:
                print(f"{FAIL} QC FAIL: {row.Sample_ID} - the following paths do not exist: {', '.join(missing_paths)}.")
                bad_sample_ids = True
    
            else:
                # Extract sequencing metadata (instrument/run info) from FASTQ header
                try:
                    result = subprocess.run(f"zcat {fq1_filepath} | head -n1", shell=True, capture_output=True, text=True, check=True)
                    fastq_header = result.stdout.strip().split("@")[1]  # e.g., "INSTRUMENT:RUN:FLOWCELL:LANE..."
                    instrument_name, run_id, flowcell_id, flowcell_lane = fastq_header.split(":")[0:4]
                    
                    # Add sample entry depending on sequencing method
                    if row['seq_method'] in ['riptide', 'LSW_lcwgs']:
                        samples_variables.append([
                            row.Sample_ID, row.seq_method, row.Library_ID, row.Sample_Barcode,
                            fq1_filepath, fq2_filepath,
                            instrument_name, run_id, flowcell_id, flowcell_lane,
                            output_dir
                        ])

                    elif row.seq_method == "gbs":
                        samples_variables.append([
                            row.Sample_ID, row.seq_method, row.Library_ID, row.Sample_Barcode,
                            fq1_filepath, "NA",  # No R2 for ddgbs
                            "NA", "NA", "GBS", "NA",
                            output_dir
                        ])

                    else:
                        print(f"{FAIL} ERROR: Invalid sequencing method {row.seq_method} for {row.Sample_ID}")
                

                except Exception as e:
                    print(f"{FAIL} QC FAIL: {row.Sample_ID} - could not extract FASTQ header info from {fq1_filepath}. {e}")
                    bad_sample_ids = True


    # Convert list to DataFrame and save as CSV
    sample_processing_input_df = pd.DataFrame(
        samples_variables,
        columns=[
            'Sample_ID', 'seq_method', 'Library_ID', 'Sample_Barcode',
            'demuxed_fq1_filepath', 'demuxed_fq2_filepath',
            'instrument_name', 'run_id', 'flowcell_id', 'flowcell_lane',
            'output_directory'
        ]
    )
    
    if bad_sample_ids: 
        print(f"{FAIL} QC FAIL: There are Sample_IDs in metadata_filepath but not in the exclude_sampleids_list, that do not have valid sample processing input and output filepaths. Will not be saving sample processing input file to {output_filepath}.")
    else:
        print(f"{PASS} QC PASS: All Sample_IDs in metadata_filepath but not in the exclude_sampleids_list, have valid sample processing input and output filepaths. Saving sample processing input file to {output_filepath}.")
        sample_processing_input_df.to_csv(output_filepath, index=True, index_label="idx")



##############################################################################################################################################################################

def create_samtools_depth_bed(positions_directory, output_directory):

    chrom_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                    "15","16","17","18","19","20","X", "Y", "M"]

    chrom_pos_df_dict = {}
    # import files manually into a dictionary
    for chrom in chrom_list:
        chrom_df = pd.read_csv(f"{positions_directory}/chr{chrom}_STITCH_baud_pos", sep="\t", 
                            header=None, names=["chrom", "pos", "ref", "alt"],
                            dtype={"pos": int})
        # chrom_df["pos_plusone"] = chrom_df["pos"] + 1
        # bed_df = chrom_df[["chrom", "pos", "pos_plusone"]]

        chrom_df["pos_minusone"] = chrom_df["pos"] - 1
        bed_df = chrom_df[["chrom", "pos_minusone", "pos"]]        

        bed_df.to_csv(f"{output_directory}/chr{chrom}_variant.bed", header=None, index=False, sep="\t")


##############################################################################################################################################################################

def create_stitch_lists(metadata_filepath, flowcells_directory, genome_shorthand, current_round, stitch_samples_filepath, output_directory, log_directory, exclude_sampleids_list, dtype_mapping):
    """
    Creates sex-specific sample lists and bam lists for STITCH runs.
    Validates the bam and bai filespath of samples to be included in these lists.
    Checks if all samples in the multi-flowcell sheet have been accounted for (sample is either excluded or include in STITCH run)
    
    Parameters:
    - metadata_filepath (str): Path to the CSV containing metadata for all samples.
    - flowcells_directory (str): Root directory where flowcell data is stored.
    - genome_shorthand (str): Identifier for the genome reference used (e.g., 'GRCr8').
    - current_round (str): Identifier for the current round (e.g. 'round11_0_HD').
    - stitch_samples_filepath (str): Path to the CSV containing samples that have passed sex-QC and are to be put through STITCH (e.g. '_stitch_samples.csv').
    - output_directory (str): Path to where the sample lists and bam lists for STITCH will be saved.
    - log_directory (str): Path to where log files containing problematic samples will be saved.
    - exclude_sampleids_list (List[str]) : List of sample IDs to exclude from processing.
    - dtype_mapping (Dict[str, type]) : Dictionary specifying the data types for reading the metadata CSV.

    If any path does not exist, an error message is printed, a log file is created, 
    and neither sample nor bam list is created. Step 5 is terminated at this point.
    """
    
    # Read in multi-flowcell sheet
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Read in samples to be processed by STITCH
    stitch_samples_df = pd.read_csv(stitch_samples_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Track problematic samples
    unaccounted_samples = []
    nobam_samples = []
    noindex_samples = []
    nobam_noindex_samples = []

    # Track all Sample IDs being checked
    checked_sampleid = 0
    
    for _, row in metadata_df.iterrows():
        
        bam_filepath = f"{flowcells_directory}/{row['flowcell_id']}/{genome_shorthand}/03_sample_processing/picard_markdup/{row['Sample_ID']}_sorted_mkDup.bam"  
        bai_filepath = f"{flowcells_directory}/{row['flowcell_id']}/{genome_shorthand}/03_sample_processing/picard_markdup/{row['Sample_ID']}_sorted_mkDup.bai"  
        
        if row['Sample_ID'] not in exclude_sampleids_list:
            
            checked_sampleid += 1
            
            # Check if samples in multi-flowcell sheet and are not in exclude list, but are missing from the STITCH samples list
            sample_info = [row['Sample_ID'], row['Library_ID'], row['Sample_Project']]
            
            if row['Sample_ID'] not in stitch_samples_df['Sample_ID'].tolist():
                unaccounted_samples.append(sample_info)
            
            # Check if bam/bai filepaths of sample exists
            no_bam = not os.path.exists(bam_filepath) or (os.path.getsize(bam_filepath) == 0)
            no_index = not os.path.exists(bai_filepath) or (os.path.getsize(bai_filepath) == 0)
                
            sample_info = [row['Sample_ID'], row['Library_ID'], row['Sample_Project']]

            if no_bam and no_index:
                print(f"nobam_noindex_samples: {sample_info}")
                nobam_noindex_samples.append(sample_info)
            elif no_bam:
                print(f"no_bam: {sample_info}")                
                nobam_samples.append(sample_info)
            elif no_index:
                print(f"no_index: {sample_info}")                
                noindex_samples.append(sample_info)    
    
    # If there are problematic samples, print them out onto the console and also log it in a text file
    for bad_sample_list, bad_sample_type, message in [(unaccounted_samples, 'unaccounted_samples', 'are missing in the STITCH samples list'),
                                                      (nobam_noindex_samples, 'nobam_noindex_samples', 'do not have bam and do not have index files'),
                                                      (nobam_samples, 'nobam_samples', 'do not have bam files'),
                                                      (noindex_samples, 'noindex_samples', 'do not have index files')]:
        if len(bad_sample_list) != 0:
            print(f"QC FAIL: There are {len(bad_sample_list)}/{checked_sampleid} samples in the multi-flowcell sheet that are not in the exclude sample ids list, but {message}.\n")
            print(f"Writing list of samples that {message} to {log_directory}/nobam_noindex_samples.txt")
        
        pd.DataFrame(bad_sample_list, columns=['Sample_ID', 'Library_ID', 'Sample_Project']).to_csv(f"{log_directory}/{bad_sample_type}.txt", index=False)                               
                    
    # If there are problematic samples, do NOT create STITCH lists
    # Force the user to address the samples
    if unaccounted_samples or nobam_noindex_samples or nobam_samples or noindex_samples:
        print("\nThere are problematic samples. No STITCH lists will be created. Exiting pipeline.")
        sys.exit(0)
    
    # If there are no problematic samples left, create STITCH lists
    else:
        
        print(f"QC SUCCESS! All multi-flowcell samples have been accounted for and relevant samples have a bam and index file. Creating STITCH lists...")
        stitch_bam_df = metadata_df.merge(stitch_samples_df[['Sample_ID', 'predicted_sex']], on='Sample_ID', how='right')
        stitch_bam_df['bamlist'] = stitch_bam_df.apply(
            lambda row: f"{flowcells_directory}/{row['flowcell_id']}/{genome_shorthand}/03_sample_processing/picard_markdup/{row['Sample_ID']}_sorted_mkDup.bam",
            axis=1
        )
        
        # Create female, male, all sexes STITCH sample lists AND bam lists
        # Note, this method should create STITCH lists in the same sample id order
        female_bamlist = stitch_bam_df[stitch_bam_df['predicted_sex'] == 'F']
        female_bam_list_filepath = f"{output_directory}/{current_round}_female_bam_list.txt"
        female_sampleid_list_filepath = f"{output_directory}/{current_round}_female_sampleid_list.txt"
        female_bamlist['bamlist'].to_csv(female_bam_list_filepath, index=False, header=False)
        female_bamlist['Sample_ID'].to_csv(female_sampleid_list_filepath, index=False, header=False)        
        print(f"Writing female bam list to {female_bam_list_filepath}.")
        print(f"Writing female sample id list to {female_sampleid_list_filepath}.")
        
        male_bamlist = stitch_bam_df[stitch_bam_df['predicted_sex'] == 'M']
        male_bam_list_filepath = f"{output_directory}/{current_round}_male_bam_list.txt"
        male_sampleid_list_filepath = f"{output_directory}/{current_round}_male_sampleid_list.txt"        
        male_bamlist['bamlist'].to_csv(male_bam_list_filepath, index=False, header=False)
        male_bamlist['Sample_ID'].to_csv(male_sampleid_list_filepath, index=False, header=False)        
        print(f"Writing male bam list to {male_bam_list_filepath}.")
        print(f"Writing male sample id list to {male_sampleid_list_filepath}.")        
        
        allsexes_bamlist = stitch_bam_df[stitch_bam_df['predicted_sex'].isin(['F', 'M'])]
        allsexes_bam_list_filepath = f"{output_directory}/{current_round}_allsexes_bam_list.txt"
        allsexes_sampleid_list_filepath = f"{output_directory}/{current_round}_allsexes_sampleid_list.txt"           
        allsexes_bamlist['bamlist'].to_csv(allsexes_bam_list_filepath, index=False, header=False)
        allsexes_bamlist['Sample_ID'].to_csv(allsexes_sampleid_list_filepath, index=False, header=False)
        print(f"Writing all sexes bam list to {allsexes_bam_list_filepath}.")
        print(f"Writing all sexes sample id list to {allsexes_sampleid_list_filepath}.")             



##############################################################################################################################################################################

def create_stitch_input(chunk_positions_directory, stitch_list_directory, positions_directory, reference_panels_directory, current_round, chunk_vcfs_directory, output_directory):
    """
    Creates a STITCH input file containing STITCH parameters that can be automated (e.g., are stable and will not be changed by user).
    Each row in this input file pertains to a STITCH chromosome chunk (job array) that is coordinated by the outputs of the `run_stitch_split_chr` function.
    Delibirates the chromosome, sex-specific method (pseudoHaploid, diploid) the STITCH chunk will be processed by.
    
    Parameters:
    - chunk_positions_directory (str): Path to the CSVs produced by the `run_stitch_split_chr` function.
    - stitch_list_directory (str): Path to the CSVs produced by the `create_stitch_lists` function.
    - positions_directory (str): Path to the STITCH positions file stated in the `pipeline_args.sh` file 
    - reference_panels_directory (str): Path to the STITCH reference panel files stated in the `pipeline_args.sh` file 
    - current_round (str): Identifier for the current round (e.g. 'round11_0_HD').
    - chunk_vcfs_directory (str): Path to where the chromosome chunk VCF files produced by STITCH will be saved.
    - output_directory (str): Path to where the STITCH input file will be saved.

    If any path does not exist, an error message is printed, a log file is created, 
    and the input file is NOT created. Step 5 is terminated at this point.
    """
    
    # List of chromosome jobs to run and the STITCH lists to be used for the chromosome run
    chrom_list = [('allsexes', '1'), ('allsexes', '2'), ('allsexes', '3'), ('allsexes', '4'), ('allsexes', '5'), ('allsexes', '6'), ('allsexes', '7'), ('allsexes', '8'), ('allsexes', '9'), ('allsexes', '10'),
              ('allsexes', '11'), ('allsexes', '12'), ('allsexes', '13'), ('allsexes', '14'), ('allsexes', '15'), ('allsexes', '16'), ('allsexes', '17'), ('allsexes', '18'), ('allsexes', '19'), ('allsexes', '20'),
              ('female', 'X'), ('male', 'X'), ('male', 'Y'), ('allsexes', 'M')]
    
    
    ########## Validate stitch input files ##########
    missing_stitch_input_files = []
    
    # Validate chromosome chunk filepaths created by the `run_stitch_split_chr` function
    for _, chrom_job in enumerate(chrom_list):
        chunk_path = f"{chunk_positions_directory}/chr{chrom_job[1]}_chunks.txt"
        if not os.path.exists(chunk_path):
            missing_stitch_input_files.append(chunk_path)
            
    
    # Validate chromosome position and reference panel files
    for _, chrom_job in enumerate(chrom_list):
        chrom_pos_filepath = f"{positions_directory}/chr{chrom_job[1]}_STITCH_baud_pos"
        chrom_refpanel_samples_filepath = f"{reference_panels_directory}/vcf1_filtered_chr{chrom_job[1]}.samples"
        chrom_refpanel_legend_filepath = f"{reference_panels_directory}/vcf1_filtered_chr{chrom_job[1]}.legend.gz"
        chrom_refpanel_hap_filepath = f"{reference_panels_directory}/vcf1_filtered_chr{chrom_job[1]}.hap.gz"
        
        pos_refpanel_filepaths = [chrom_pos_filepath, chrom_refpanel_samples_filepath, chrom_refpanel_legend_filepath, chrom_refpanel_hap_filepath]
        
        for pos_refpanel_filepath in pos_refpanel_filepaths:
            if not os.path.exists(pos_refpanel_filepath) or (os.path.getsize(pos_refpanel_filepath) == 0):
                missing_stitch_input_files.append(pos_refpanel_filepath)        
    
    # Validate stitch list files
    stitch_list_file_suffixes = ["allsexes_bam_list.txt", "allsexes_sampleid_list.txt",
                                 "male_bam_list.txt", "male_sampleid_list.txt",
                                 "female_bam_list.txt", "female_sampleid_list.txt"]
    
    for stitch_list_file_suffix in stitch_list_file_suffixes:
        stitch_list_filepath = f"{stitch_list_directory}/{current_round}_{stitch_list_file_suffix}"
        if not os.path.exists(stitch_list_filepath) or (os.path.getsize(stitch_list_filepath) == 0):
            missing_stitch_input_files.append(stitch_list_filepath)
            
    # Validate STITCH output directory
    if not os.path.exists(chunk_vcfs_directory):
        missing_stitch_input_files.append(chunk_vcfs_directory)
  
    # If there are missing STITCH list files, do NOT create STITCH input file
    # Force the user to address the missing files            
    if missing_stitch_input_files:
        print("QC FAIL: The following STITCH input files are missing:")
        for missing_file in missing_stitch_input_files:
            print(missing_file)
            
        print("No STITCH input file will be created. Exiting pipeline.")
        sys.exit(1)
    
    else:
        print("QC SUCCESS! All STITCH input filepaths have been validated. Creating STITCH input file...")
    
    ########## Create tuples of the starting and ending position of each STITCH chunk ##########
    allchrom_tuple_dict = {}
    
    # Iterate through the list of chromosome jobs and create chromosome chunks for each
    for _, chrom_job in enumerate(chrom_list):
        chunk_df = pd.read_csv(f"{chunk_positions_directory}/chr{chrom_job[1]}_chunks.txt", header=None, names=['chunk_positions'])

        # Store the chromosome chunks/tuples
        chrom_tuple_list = []

        for index, row in chunk_df.iterrows():

            # Defines a chunk, as starting at the position stated by the current row, and ending at the position stated by the next row
            # The number of chunks of a chromosome is then n - 1, where n is the number of rows in the `_chunks.txt`
            # The tuple_val variable is a tuple (chromosome, starting position of chunk, ending position of chunk)
            if index != (chunk_df.shape[0] - 1):
                tuple_val = (f"chr{chrom_job[1]}", row['chunk_positions'], chunk_df.loc[index+1, 'chunk_positions'])
                chrom_tuple_list.append(tuple_val)
        
        # Append the list of tuples (for the chromosome job) into a dictionary using a chromosome-key
        allchrom_tuple_dict[f"chr{chrom_job[1]}"] = chrom_tuple_list
         
    
    ########## Create a dataframe coordinating STITCH chunks ##########
    
    # Store dataframes for each chromosome job
    allchrom_df_list = []
    
    # Iterate through the list of chromosome jobs and create chromosome chunks for each
    for _, chrom_job in enumerate(chrom_list):
        
        # Create a dataframe from the list of tuples of a chromosome job
        # Each row pertains to a chromosome chunk, and each column is a stable STITCH parameter
        chrom_df = pd.DataFrame(allchrom_tuple_dict[f"chr{chrom_job[1]}"], columns=['chr', 'regionStart', 'regionEnd'])
        chrom_df['sex'] = chrom_job[0]
        chrom_df['k'] = 8
        chrom_df['nGen'] = 100
        chrom_df['niterations'] = 2
        chrom_df['posfile'] = f"{positions_directory}/chr{chrom_job[1]}_STITCH_baud_pos"
        chrom_df['reference_panels'] = f"{reference_panels_directory}/vcf1_filtered_chr{chrom_job[1]}"

        # State the method type, bam list filepath, and sample id list filepath depending on sex and chromosome
        if chrom_job[0] == 'allsexes':
            chrom_df['method'] = 'diploid'
            chrom_df['bamlist'] = f"{stitch_list_directory}/{current_round}_allsexes_bam_list.txt"
            chrom_df['sampleNames_file'] = f"{stitch_list_directory}/{current_round}_allsexes_sampleid_list.txt"     
        
        elif chrom_job[0] == 'male':
            # chrom_df['method'] = 'pseudoHaploid'
            chrom_df['method'] = 'diploid-inbred'            
            chrom_df['bamlist'] = f"{stitch_list_directory}/{current_round}_male_bam_list.txt"
            chrom_df['sampleNames_file'] = f"{stitch_list_directory}/{current_round}_male_sampleid_list.txt"   

        elif chrom_job[0] == 'female':
            chrom_df['method'] = 'diploid'  
            chrom_df['bamlist'] = f"{stitch_list_directory}/{current_round}_female_bam_list.txt"
            chrom_df['sampleNames_file'] = f"{stitch_list_directory}/{current_round}_female_sampleid_list.txt"
            
        if chrom_job == ('allsexes', 'M'):
            # chrom_df['method'] = 'pseudoHaploid'  
            chrom_df['method'] = 'diploid-inbred'              
            chrom_df['bamlist'] = f"{stitch_list_directory}/{current_round}_allsexes_bam_list.txt"
            chrom_df['sampleNames_file'] = f"{stitch_list_directory}/{current_round}_allsexes_sampleid_list.txt"        
        
        # State the output filepath of the STITCH chromosome chunk run (i.e., job array)
        # Create the output directory of the STITCH chromosome chunk run (i.e., job array)
        chrom_df['stitch_path'] = None
        for index, row in chrom_df.iterrows():
            regionStart_filename = row['regionStart'] + 1
            chrom_df.loc[index, 'output_filepath'] = f"{chunk_vcfs_directory}/{row['chr']}_{row['sex']}/stitch.{row['chr']}.{regionStart_filename}.{row['regionEnd']}.vcf.gz"
            chrom_df.loc[index, 'output_directory'] = f"{chunk_vcfs_directory}/{row['chr']}_{row['sex']}"
            # create output directory
            os.makedirs(f"{chunk_vcfs_directory}/{row['chr']}_{row['sex']}", exist_ok=True)
            
        # Rearrange the columns such that the parameters are mostly in order of the `STITCH.R` script
        chrom_df = chrom_df[['chr', 'sex', 'k', 'nGen', 'niterations', 'method',
                             'output_filepath', 'output_directory',                          
                             'bamlist', 'sampleNames_file',
                             'posfile', 'regionStart', 'regionEnd',
                             'reference_panels'
                            ]]
        
        # Create an 'idx' column to be referenced by the `create_stitch_arrays` function
        # The `idx` values dictate the array values to be run
        # chrom_df = chrom_df.reset_index(names=['idx'])
        allchrom_df_list.append(chrom_df)
        
    # Concatenate all the dataframes for each chromosome job into a master stitch input dataframe
    stitch_input = pd.concat(allchrom_df_list).reset_index(drop=True)
    stitch_input = stitch_input.reset_index(names=['idx'])
        
    # Create a summary dataframe of the stitch input for the user perusal
    stitch_input_summary = stitch_input[['chr', 'sex', 'method', 'bamlist', 'sampleNames_file']].value_counts(dropna=False).reset_index()
    
    # Export the stitch input and summary dataframes
    stitch_input_filepath = f"{output_directory}/{current_round}_stitch_input.csv"
    stitch_input_summary_filepath = f"{output_directory}/{current_round}_stitch_input_summary.csv"
    stitch_input.to_csv(stitch_input_filepath, index=False)
    stitch_input_summary.to_csv(stitch_input_summary_filepath, index=False)
    print(f"Writing STITCH input file to {stitch_input_filepath}.")
    print(f"Writing STITCH input summary file to {stitch_input_summary_filepath}.")          
    
    
    return stitch_input, stitch_input_summary
    

##############################################################################################################################################################################

def create_pedigree_file(metadata_filepath, imputed_sex_filepath, current_round, exclude_sampleids_list, output_directory, dtype_mapping):
    
    # Read in metadata_filepath
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    metadata_df = metadata_df.drop(columns=['sex'])

    # Read in imputed sex
    imputed_sex_df = pd.read_csv(imputed_sex_filepath, dtype=dtype_mapping)[['Sample_ID', 'predicted_sex']]
    imputed_sex_df = imputed_sex_df.rename(columns={'predicted_sex' : 'sex'})
    
    # Remove bad samples
    metadata_df = metadata_df[~metadata_df['Sample_ID'].isin(exclude_sampleids_list)]

    # Create a DataFrame mapping (potential) dam names to their sample IDs
    dam_df = metadata_df[['Sample_Name', 'Sample_ID']].rename(columns={'Sample_Name': 'dams', 'Sample_ID': 'dams_sample_id'})

    # Create a DataFrame mapping (potential) sire names to their sample IDs
    sire_df = metadata_df[['Sample_Name', 'Sample_ID']].rename(columns={'Sample_Name': 'sires', 'Sample_ID': 'sires_sample_id'})

    # Merge dam and sire sample IDs back into the original DataFrame based on 'dams' and 'sires' names
    parental_sampleid_merge = metadata_df.merge(dam_df, on='dams', how='left').merge(sire_df, on='sires', how='left')

    # Merge sex data into the dataframe
    sex_merged_df = parental_sampleid_merge.merge(imputed_sex_df, on='Sample_ID', how='right')

    # Create pedigree dataframe (select relevant columns)
    pedigree = sex_merged_df[['Sample_ID', 'dams_sample_id', 'sires_sample_id', 'sex']]

    # Rename pedigree columns to be MER script compatible
    # The MER script will use the columns 'id', 'mother', 'father'
    pedigree_renamed = pedigree.rename(columns={'Sample_ID':'id', 'dams_sample_id': 'mother', 'sires_sample_id': 'father'})

    # Drop rows where either mother or father data is missing
    pedigree_renamed = pedigree_renamed.dropna(axis=0, subset=['mother', 'father'])

    pedigree_renamed.to_csv(f"{output_directory}/{current_round}_pedigree.csv", header=True, index=False)

