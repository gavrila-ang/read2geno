 #!/usr/bin/env python
import pandas as pd
import os
import csv
import warnings
import subprocess
import sys
from tqdm import tqdm
from utils.validate_paths import validate_paths
from utils.print_aesthetics import *
warnings.simplefilter(action='ignore', category=FutureWarning)


##############################################################################################################################################################################


def create_demux_arrays(demux_input_filepath, metadata_filepath, flowcells_directory, genome_shorthand, exclude_sampleids_list, summary_filepath, output_filepath, dtype_mapping):
    """
    Identifies libraries containing samples that failed demultiplexing and generates
    an array of their indices from the demux input file for reprocessing.

    Parameters:
    - demux_input_filepath (str): Path to the CSV file with demux input data.
    - metadata_filepath (str): Path to the flowcell-level sample sheet.
    - exclude_sampleids_list (list): Sample IDs to ignore when checking for failures.
    - output_filepath (str): Path to write the output array of failed library indices.
    - dtype_mapping (dict): Dictionary specifying dtypes for reading CSVs.
    
    Returns:
    - int: Total count of undemultiplexed Sample_IDs across all libraries.
    """

    # Validate that required files and directories exist
    output_directory = os.path.dirname(output_filepath)
    validate_paths([demux_input_filepath, metadata_filepath, output_directory])
    
    # Load input data
    demux_input = pd.read_csv(demux_input_filepath, dtype=dtype_mapping).reset_index(drop=True)
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Track total number of undemuxed samples
    n_undemuxed_sampleids = 0
    
    # Track total number of samples to be demuxed
    n_total_sampleids_to_demux = 0
    
    # Dictionary to track undemuxed sample IDs per library
    undemuxed_sampleids_dict = {lib: [] for lib in metadata_df['Library_ID'].unique()}
    undemuxed_sampleids_list = []

    # Iterate through sample sheet
    for _, row in metadata_df.iterrows():
        if row.Sample_ID not in exclude_sampleids_list:

            n_total_sampleids_to_demux += 1

            base_path = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/02_demux/demux/{row.Library_ID}"

            # Check existence of required FASTQs based on sequencing method
            if row['seq_method'] in ['riptide', 'LSW_lcwgs']:
                fq1_filepath = f"{base_path}/{row.Sample_ID}-{row.Sample_Name}-{row.Sample_Barcode}_R1.fastq.gz"
                fq2_filepath = f"{base_path}/{row.Sample_ID}-{row.Sample_Name}-{row.Sample_Barcode}_R2.fastq.gz"
                if not os.path.exists(fq1_filepath) or not os.path.exists(fq2_filepath) \
                    or os.path.getsize(fq1_filepath) == 0 or os.path.getsize(fq2_filepath) == 0:
                    undemuxed_sampleids_dict[row.Library_ID].append(row.Sample_ID)
                    undemuxed_sampleids_list.append(row.Sample_ID)
                    n_undemuxed_sampleids += 1

            elif row.seq_method == "gbs":
                fq1_filepath = f"{base_path}/{row.Sample_ID}.fq.gz"
                if not os.path.exists(fq1_filepath) or os.path.getsize(fq1_filepath) == 0:
                    undemuxed_sampleids_dict[row.Library_ID].append(row.Sample_ID)
                    undemuxed_sampleids_list.append(row.Sample_ID)                    
                    n_undemuxed_sampleids += 1

            else:
                print(f"ERROR: Invalid sequencing method {row.seq_method} for {row.Sample_ID}")

    # Determine which libraries had undemuxed samples and collect their indices
    rerun_arrays = []
    for library_id, sample_list in undemuxed_sampleids_dict.items():
        
        flowcell_id = metadata_df.query(f'Library_ID == "{library_id}"')["flowcell_id"].unique().item()
        
        if library_id in demux_input.Library_ID.tolist(): # excludes libraries in flowcell sheet but omitted from demux_input         

            library_idx = demux_input.loc[demux_input.Library_ID == library_id, "idx"].item()

            if sample_list:  # Non-empty list means failures exist
                print(f"{FAIL} QC FAIL: Flowcell ID {flowcell_id}, Library {library_id} has {len(sample_list)} Sample_IDs not demuxed (index {library_idx}).")
                rerun_arrays.append(library_idx)
            else:
                print(f"{PASS} QC PASS: Flowcell ID {flowcell_id}, Library {library_id} has all Sample_IDs demuxed (index {library_idx}).")

        else:
            print(f"{WARN} QC WARNING: Flowcell ID {flowcell_id}, Library {library_id} was not included in the demux_input file due to missing filepaths.")


    undemuxed_sampleids_df = metadata_df.loc[metadata_df['Sample_ID'].isin(undemuxed_sampleids_list), ['Sample_ID','Library_ID', 'flowcell_id', 'Sample_Project']]
    print(f"Saving undemuxed sample ids file to {summary_filepath}.")
    undemuxed_sampleids_df.to_csv(summary_filepath, index=None)

    # Write the rerun array indices to a CSV
    print(f"Saving demuxing array file to {output_filepath}.")
    with open(output_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for idx in rerun_arrays:
            writer.writerow([idx])

    # Return total number of undemuxed samples
    return n_undemuxed_sampleids, n_total_sampleids_to_demux
    
    

##############################################################################################################################################################################

def create_sample_processing_arrays(metadata_filepath, sample_processing_input_filepath, flowcells_directory, genome_shorthand, exclude_sampleids_list,
                                    output_filepath, dtype_mapping):
    """
    Identifies samples with incomplete processing by checking the presence of expected output files
    for each step in the pipeline. Generates an array indicating which samples need to be rerun and 
    from which step.

    Parameters:
    - sample_processing_input_filepath (str): Path to the sample processing input CSV file.
    - exclude_sampleids_list (list): List of Sample_IDs to skip from evaluation.
    - output_filepath (str): Path to write the array of samples to be rerun.
    - dtype_mapping (dict): Dictionary specifying data types for reading CSV.

    Returns:
    - int: Total number of unprocessed or partially processed samples.
    """

    # Validate paths
    output_directory = os.path.dirname(output_filepath)
    validate_paths([sample_processing_input_filepath, output_directory])
    
    # Read sample processing input
    sample_processing_input_df = pd.read_csv(sample_processing_input_filepath, dtype=dtype_mapping).reset_index(drop=True)
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # List to collect samples that need to be rerun
    rerun_arrays = []  

    # Track total number of samples that are not fully processed
    n_unprocessed_sampleids = 0 

    # Track total number of samples to be processed
    n_total_sampleids_to_process = 0

    # Iterate through sample sheet
    for _, row in metadata_df.iterrows():

        if row.Sample_ID not in exclude_sampleids_list:

            n_total_sampleids_to_process += 1
            
            base_path = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing"

            # Define expected output files per step
            # Note the samtools depth chrM files are not checked because they could naturally be empty
            output_filepaths_dict = {
                "bbduk": [
                    f"{base_path}/bbduk/{row.Sample_ID}_adapter_trimmed_R1.fastq.gz",
                    f"{base_path}/bbduk/{row.Sample_ID}_adapter_trimmed_R2.fastq.gz"
                ],
                "cutadapt": [
                    f"{base_path}/cutadapt/{row.Sample_ID}_qualen_trimmed_R1.fastq.gz",
                    f"{base_path}/cutadapt/{row.Sample_ID}_qualen_trimmed_R2.fastq.gz"
                ],
                "fastqc_paired": [
                    f"{base_path}/fastqc_trim/{row.Sample_ID}_adapter_trimmed_R1_fastqc.zip",
                    f"{base_path}/fastqc_trim/{row.Sample_ID}_adapter_trimmed_R2_fastqc.zip",
                    f"{base_path}/fastqc_trim/{row.Sample_ID}_qualen_trimmed_R1_fastqc.zip",
                    f"{base_path}/fastqc_trim/{row.Sample_ID}_qualen_trimmed_R2_fastqc.zip"
                ],
                "fastqc_single": [f"{base_path}/fastqc_trim/{row.Sample_ID}_fastqc.zip"],              
                "bwamem": [f"{base_path}/bwamem/{row.Sample_ID}.sam"],
                "samtools_sort": [f"{base_path}/samtools_sort/{row.Sample_ID}_sorted.bam"],
                "picard_markdup": [
                    f"{base_path}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam", 
                    f"{base_path}/picard_markdup/{row.Sample_ID}_sorted_mkDup.metrics"
                ],
                "samtools_index": [f"{base_path}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bai"],
                "samtools_flagstat": [f"{base_path}/samtools_flagstat/{row.Sample_ID}.flagstat"],
                "samtools_idxstats": [f"{base_path}/samtools_idxstats/{row.Sample_ID}.idxstats"],
                "mosdepth": [f"{base_path}/mosdepth/{row.Sample_ID}.mosdepth.global.dist.txt", 
                             f"{base_path}/mosdepth/{row.Sample_ID}.mosdepth.summary.txt"],
                "mosdepth_exclude_duplicates": [f"{base_path}/mosdepth/{row.Sample_ID}_excluded_duplicates.mosdepth.global.dist.txt", 
                                                f"{base_path}/mosdepth/{row.Sample_ID}_excluded_duplicates.mosdepth.summary.txt"],
                "mosdepth_exclude_unmapped": [f"{base_path}/mosdepth/{row.Sample_ID}_excluded_unmapped.mosdepth.global.dist.txt", 
                                                f"{base_path}/mosdepth/{row.Sample_ID}_excluded_unmapped.mosdepth.summary.txt"],                                                
                "samtools_depth": [
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr1.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr2.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr3.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr4.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr5.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr6.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr7.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr8.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr9.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr10.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr11.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr12.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr13.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr14.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr15.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr16.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr17.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr18.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr19.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chr20.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chrX.depth",
                    f"{base_path}/samtools_depth/{row.Sample_ID}_chrY.depth"
                    # f"{base_path}/samtools_depth/{row.Sample_ID}_chrM.depth"
                ],
                "samtools_stats": [f"{base_path}/samtools_stats/{row.Sample_ID}.stats"],
                "samtools_view": [
                    f"{base_path}/samtools_view/{row.Sample_ID}_chrX_count.csv",
                    f"{base_path}/samtools_view/{row.Sample_ID}_chrY_count.csv",
                    f"{base_path}/samtools_view/{row.Sample_ID}_allchrom_count.csv"
                ],
                "cutadapt_5": [f"{base_path}/cutadapt/{row.Sample_ID}_adapter_trimmed.fq.gz"],
                "cutadapt_3": [f"{base_path}/cutadapt/{row.Sample_ID}_qualen_trimmed.fq.gz"]
            }

            # Define the pipeline steps per sequencing method
            substeps = ["picard_markdup", "samtools_index", "samtools_flagstat", "samtools_view", "mosdepth_exclude_duplicates"]

            # for substep in substeps:
            #     outfiles = output_filepaths_dict[substep]
            #     for outfile in outfiles:
            #         if not os.path.exists(outfile) and not os.path.getsize(outfile) > 0:
            #             print(outfile)
                        
            # Enumerate the substeps starting from 1 (e.g., 1, 2, 3...) and convert to a list
            # This keeps substep numbering consistent with external references
            enumerated_substeps = list(enumerate(substeps, start=1))

            restart_from_substep = None
            failed_substep = None

            # Iterate over the substeps to find the *last successfully completed* substep
            for step_index, step_name in enumerated_substeps:

                # Get the list of expected output files for the current substep
                expected_files = output_filepaths_dict.get(step_name, [])
                
                # If all expected files for this substep exist and are non-empty, assume it completed successfully
                if all(os.path.exists(path) and os.path.getsize(path) > 0 for path in expected_files):
                    continue
                
                # If not, the current substep (unless first substep) 
                # likely failed because of a truncated output from the previous substep
                # and thus mark the previous substep for rerun
                else:
                    # restart_from_substep = step_index - 1 if step_index != 1 else step_index
                    failed_substep = step_name
                    restart_from_substep = 1 # uncomment if things aren't going well 
                    break

            # # If any substep has failed, mark the sample for rerun
            if restart_from_substep is not None:
                # print(f"QC FAIL: Sample_ID '{row.Sample_ID}' from Library_ID '{row.Library_ID}' is restarting from step 3.{restart_from_substep}")
                # Look up the index ('idx') of the current sample in the sample_processing_input_df
                # This index corresponds to the sample's position in the input array for rerun tracking
                try:
                    spi_idx = sample_processing_input_df.loc[sample_processing_input_df.Sample_ID == row.Sample_ID, "idx"].item()
                except Exception as e:
                    print(row.Sample_ID, e)
                # Append the sample's rerun information:
                rerun_arrays.append([spi_idx, row.Library_ID, row.Sample_ID, failed_substep, restart_from_substep])
                # Increment the count of samples needing reprocessing
                n_unprocessed_sampleids += 1

    # Write rerun array to CSV
    print(f"Saving sample processing array file to {output_filepath}.")
    with open(output_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["idx", "Library_ID", "Sample_ID", "failed_substep_name", "restart_from_substep"])
        writer.writerows(rerun_arrays)

    return n_unprocessed_sampleids, n_total_sampleids_to_process



##############################################################################################################################################################################

def create_stitch_arrays(metadata_filepath, round_directory, current_round, stitch_input_filepath, output_filepath, exclude_sampleids_list, dtype_mapping):
    """
    Creates a file containing three columns [idx, chrom, sex].
    The `chrom` and `sex` columns contains values stating the chromosome JOBS to be run. The number of chromosome JOBs to be run is stated by the chom_list.
    The `idx` column contains values ranging from 0 to number of chromosome chunks for each chromosome. In other words, each `idx` value is a chromosome chunk, aka, an array job.
    Each chromosome JOB will access rows pertaining to that job using the `chrom` column.
    The ARRAYS to be run/rerun, for a given chromosome JOB, will be dictated by the `idx` column.
    
    Parameters:
    - metadata_filepath (str): Path to the CSV containing metadata for all samples.
    - stitch_input_filepath (str): Path to the CSVs produced by the `create_stitch_input` function.
    - output_filepath (str): Path to where the STITCH array file will be saved.    
    - exclude_sampleids_list (List[str]) : List of sample IDs to exclude from processing.    
    - dtype_mapping (Dict[str, type]) : Dictionary specifying the data types for reading the metadata CSV.
    """
    
    # Read in the stitch_input file
    stitch_input_df = pd.read_csv(stitch_input_filepath)
    n_tostitch_chunks = stitch_input_df.shape[0]

    # Read in multi-flowcell sheet
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # Read in sample id lists
    allsexes_sid = pd.read_csv(f"{round_directory}/output/05_stitch/stitch_lists/{current_round}_allsexes_sampleid_list.txt", header=None)
    female_sid = pd.read_csv(f"{round_directory}/output/05_stitch/stitch_lists/{current_round}_female_sampleid_list.txt", header=None)
    male_sid = pd.read_csv(f"{round_directory}/output/05_stitch/stitch_lists/{current_round}_male_sampleid_list.txt", header=None)
    
    # List of chromosome jobs to run and the STITCH lists to be used for the chromosome run
    # Iterate through this list rather than row-by-row of stitch input df to ensure all chromosomes are processed
    chrom_list = [('allsexes', '1'), ('allsexes', '2'), ('allsexes', '3'), ('allsexes', '4'), ('allsexes', '5'), ('allsexes', '6'), ('allsexes', '7'), ('allsexes', '8'), ('allsexes', '9'), ('allsexes', '10'),
              ('allsexes', '11'), ('allsexes', '12'), ('allsexes', '13'), ('allsexes', '14'), ('allsexes', '15'), ('allsexes', '16'), ('allsexes', '17'), ('allsexes', '18'), ('allsexes', '19'), ('allsexes', '20'),
              ('female', 'X'), ('male', 'X'), ('male', 'Y'), ('allsexes', 'M')]
    
    # Track the chromosome chunks/arrays to be rerun
    rerun_array_list = []
    
    # Count how many chromosome chunks remain unprocessed
    n_unstitched_allchrom_chunks = 0 
    
    # Iterate through the list of chromosome jobs and 
    for _, chrom_job in enumerate(chrom_list):
        
        n_unstitched_chrom_chunks = 0
        chrom_stitch_input_df = stitch_input_df[(stitch_input_df['chr'] == f"chr{chrom_job[1]}") & (stitch_input_df['sex'] == chrom_job[0])]
        
        for _, row in chrom_stitch_input_df.iterrows():
            
            # Determine if a VCF file for the chromosome chunk has been produced by STITCH; if not output the corresponding idx/array values.
            if not os.path.exists(row['output_filepath']) or (os.path.getsize(row['output_filepath']) == 0):
                rerun_tuple = (row['idx'], chrom_job[1], chrom_job[0])
                rerun_array_list.append(rerun_tuple)
                n_unstitched_allchrom_chunks += 1
                n_unstitched_chrom_chunks += 1

            elif not os.path.exists(f"{row['output_filepath']}.tbi") or (os.path.getsize(f"{row['output_filepath']}.tbi") == 0):
                rerun_tuple = (row['idx'], chrom_job[1], chrom_job[0])
                rerun_array_list.append(rerun_tuple)
                n_unstitched_allchrom_chunks += 1
                n_unstitched_chrom_chunks += 1       

            else:
                # If a VCF file has been produced, determine if all the relevant samples are present in the VCF file
                # Relevant samples defined as sampels in the multi-flowcell sheet but not in the exclude samples list
                rerun_tuple = (row['idx'], chrom_job[1], chrom_job[0])
                if row['sex'] == 'allsexes':
                    chunk_file_samples = subprocess.run(f"zcat {row['output_filepath']} | grep -m1 '#CHROM'", shell=True, text=True, capture_output=True).stdout.strip().split('\t')[9:]
                    # relevant_samples = set(metadata_df['Sample_ID'].tolist()).difference(set(exclude_sampleids_list))
                    relevant_samples = set(allsexes_sid[0].tolist()).difference(set(exclude_sampleids_list))
                    sample_missing_in_chunk = relevant_samples.difference(set(chunk_file_samples))
                    total_samples_in_chunk = len(relevant_samples)
                elif row['sex'] == 'female':
                    chunk_file_samples = subprocess.run(f"zcat {row['output_filepath']} | grep -m1 '#CHROM'", shell=True, text=True, capture_output=True).stdout.strip().split('\t')[9:]
                    # female_sampleids = metadata_df.loc[metadata_df['sex'] == 'F', 'Sample_ID'].tolist()
                    # relevant_samples = set(female_sampleids).difference(set(exclude_sampleids_list))
                    relevant_samples = set(female_sid[0].tolist()).difference(set(exclude_sampleids_list))
                    sample_missing_in_chunk = relevant_samples.difference(set(chunk_file_samples))      
                    total_samples_in_chunk = len(relevant_samples)
                elif row['sex'] == 'male':
                    chunk_file_samples = subprocess.run(f"zcat {row['output_filepath']} | grep -m1 '#CHROM'", shell=True, text=True, capture_output=True).stdout.strip().split('\t')[9:]
                    # male_sampleids = metadata_df.loc[metadata_df['sex'] == 'M', 'Sample_ID'].tolist()
                    # relevant_samples = set(male_sampleids).difference(set(exclude_sampleids_list))
                    relevant_samples = set(male_sid[0].tolist()).difference(set(exclude_sampleids_list))
                    sample_missing_in_chunk = relevant_samples.difference(set(chunk_file_samples))      
                    total_samples_in_chunk = len(relevant_samples)                    
                    
                if sample_missing_in_chunk:
                    rerun_array_list.append(rerun_tuple)
                    n_unstitched_allchrom_chunks += 1
                    n_unstitched_chrom_chunks += 1
                    print(f"For chromosome job {row['chr']}, {row['sex']} there are {len(sample_missing_in_chunk)}/{total_samples_in_chunk} relevant samples missing in the STITCH file {row['output_filepath']}.")
    
        if n_unstitched_chrom_chunks != 0:
            print(f"For chromosome job {row['chr']}, {row['sex']} there are {n_unstitched_chrom_chunks}/{chrom_stitch_input_df.shape[0]} chromosome chunks to be run.")

    # If there are arrays left to be rerun, export the rerun array dataframe
    if n_unstitched_allchrom_chunks != 0 :
        print(f"\nQC FAIL: There are {n_unstitched_allchrom_chunks}/{n_tostitch_chunks} chunks in total to run.")
        print(f"Writing STITCH array file to {output_filepath}.")          
    else:
        print(f"\nQC SUCCESS! There are no chunks left in total to run.")

    rerun_df = pd.DataFrame(rerun_array_list, columns=['idx', 'chr', 'sex'])
    rerun_df.to_csv(output_filepath, index=False)

    return n_unstitched_allchrom_chunks, n_tostitch_chunks