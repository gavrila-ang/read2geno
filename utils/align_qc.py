#!/usr/bin/env python
import pandas as pd
import os
import subprocess
import numpy as np
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



##############################################################################################################################################################################

def bam_bai_check(metadata_filepath, exclude_sampleids_list, genome_shorthand, flowcells_directory, log_directory, dtype_mapping):

    # Read in multi-flowcell sheet
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Track all Sample IDs being checked
    checked_sampleid = 0
    
    # Track which Sample IDs don't have a bam and bai file required for STITCH
    nobam_samples = []
    noindex_samples = []
    nobam_noindex_samples = []
    
    for _, row in metadata_df.iterrows():
        
        bam_filepath = f"{flowcells_directory}/{row['flowcell_id']}/{genome_shorthand}/03_sample_processing/picard_markdup/{row['Sample_ID']}_sorted_mkDup.bam"  
        bai_filepath = f"{flowcells_directory}/{row['flowcell_id']}/{genome_shorthand}/03_sample_processing/picard_markdup/{row['Sample_ID']}_sorted_mkDup.bai"  
        
        if row['Sample_ID'] not in exclude_sampleids_list:
            
            checked_sampleid += 1
            
            # Check if bam/bai filepaths of sample exists
            no_bam = not os.path.exists(bam_filepath) or (os.path.getsize(bam_filepath) == 0)
            no_index = not os.path.exists(bai_filepath) or (os.path.getsize(bai_filepath) == 0)
                
            sample_info = [row['Sample_ID'], row['Library_ID'], row['Sample_Project']]

            if no_bam and no_index:
                nobam_noindex_samples.append(sample_info)
            elif no_bam:
                nobam_samples.append(sample_info)
            elif no_index:
                noindex_samples.append(sample_info)    

    if nobam_noindex_samples:
        print(f"There are {len(nobam_noindex_samples)}/{checked_sampleid} samples in the multi-flowcell sheet that are not in the exclude sample ids list, with NO bam and NO index files.")
        print(f"Writing list of samples with no bam and no index files to {log_directory}/nobam_noindex_samples.txt")
    pd.DataFrame(nobam_noindex_samples, columns=['Sample_ID', 'Library_ID', 'Sample_Project']).to_csv(f"{log_directory}/nobam_noindex_samples.txt", index=False)
        
    if noindex_samples:
        print(f"There are {len(noindex_samples)}/{checked_sampleid} samples in the multi-flowcell sheet that are not in the exclude sample ids list, with NO index files.")
        print(f"Writing list of samples with no index files to {log_directory}/noindex_samples.txt")
    pd.DataFrame(noindex_samples, columns=['Sample_ID', 'Library_ID', 'Sample_Project']).to_csv(f"{log_directory}/noindex_samples.txt", index=False)
        
    if nobam_samples:
        print(f"There are {len(nobam_samples)}/{checked_sampleid} samples in the multi-flowcell sheet that are not in the exclude sample ids list, with NO bam files.")
        print(f"Writing list of samples with no bam files to {log_directory}/nobam_samples.txt")        
    pd.DataFrame(nobam_samples, columns=['Sample_ID', 'Library_ID', 'Sample_Project']).to_csv(f"{log_directory}/nobam_samples.txt", index=False)
        
    if not(nobam_noindex_samples or noindex_samples or nobam_samples):
        print(f"QC SUCCESS! All samples in the multi-flowcell sheet that are not in the exclude sample list have bam/bai files.")
    else:
        print(f"QC FAIL: There are samples missing bam/bai files. Exiting pipeline.")
        sys.exit(0)


##############################################################################################################################################################################

def alignment_qc_data_check(metadata_filepath, exclude_sampleids_list, genome_shorthand, flowcells_directory, dtype_mapping):
    
    """
    Checks for the existence and completeness of alignment QC files for each sample.

    Parameters
    - metadata_filepath (str) : Path to the CSV file containing metadata for all samples across flowcells. 
    - exclude_sampleids_list (list[str]) : List of Sample_IDs to be excluded from QC checks.
    - flowcells_directory (str) : Root directory where flowcell-specific subdirectories are stored.
    - genome_shorthand (str) : Identifier used to locate genome-specific QC paths (e.g., 'GRCr8').
    - dtype_mapping (dict[str, type]) : Mapping of column names to data types used when reading the CSV file.
    """
    
    # Load the multiflowcell CSV into a DataFrame
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Initialize a set to collect Sample_IDs that fail QC
    bad_samples = set()
    
    # Counter to track how many samples are checked (excluding ones in exclude sample list)
    checked_sampleid = 0

    # Iterate through each row (sample) in the DataFrame
    for _, row in metadata_df.iterrows():
        
        # Define expected output files for align qc steps
        output_files_dict = {
        "fastqc": [f"fastqc_trim/{row.Sample_ID}_adapter_trimmed_R1_fastqc.zip",
                   f"fastqc_trim/{row.Sample_ID}_adapter_trimmed_R2_fastqc.zip",
                   f"fastqc_trim/{row.Sample_ID}_qualen_trimmed_R1_fastqc.zip",
                   f"fastqc_trim/{row.Sample_ID}_qualen_trimmed_R2_fastqc.zip"],
        "samtools_flagstat": [f"samtools_flagstat/{row.Sample_ID}.flagstat"],
        "samtools_idxstats": [f"samtools_idxstats/{row.Sample_ID}.idxstats"],
        "mosdepth": [f"mosdepth/{row.Sample_ID}.mosdepth.global.dist.txt", 
                     f"mosdepth/{row.Sample_ID}.mosdepth.summary.txt"],
        "samtools_depth": [
            f"samtools_depth/{row.Sample_ID}_chr1.depth",
            f"samtools_depth/{row.Sample_ID}_chr2.depth",
            f"samtools_depth/{row.Sample_ID}_chr3.depth",
            f"samtools_depth/{row.Sample_ID}_chr4.depth",
            f"samtools_depth/{row.Sample_ID}_chr5.depth",
            f"samtools_depth/{row.Sample_ID}_chr6.depth",
            f"samtools_depth/{row.Sample_ID}_chr7.depth",
            f"samtools_depth/{row.Sample_ID}_chr8.depth",
            f"samtools_depth/{row.Sample_ID}_chr9.depth",
            f"samtools_depth/{row.Sample_ID}_chr10.depth",
            f"samtools_depth/{row.Sample_ID}_chr11.depth",
            f"samtools_depth/{row.Sample_ID}_chr12.depth",
            f"samtools_depth/{row.Sample_ID}_chr13.depth",
            f"samtools_depth/{row.Sample_ID}_chr14.depth",
            f"samtools_depth/{row.Sample_ID}_chr15.depth",
            f"samtools_depth/{row.Sample_ID}_chr16.depth",
            f"samtools_depth/{row.Sample_ID}_chr17.depth",
            f"samtools_depth/{row.Sample_ID}_chr18.depth",
            f"samtools_depth/{row.Sample_ID}_chr19.depth",
            f"samtools_depth/{row.Sample_ID}_chr20.depth",
            f"samtools_depth/{row.Sample_ID}_chrX.depth",
            f"samtools_depth/{row.Sample_ID}_chrY.depth",
            f"samtools_depth/{row.Sample_ID}_chrM.depth"
        ],
        "samtools_stats": [f"samtools_stats/{row.Sample_ID}.stats"]
        }
        
        # Skip samples that are excluded or not sequenced using 'riptide'
        if row.Sample_ID in exclude_sampleids_list or row.seq_method != "riptide":
            continue
        
        # Count this sample as checked
        checked_sampleid += 1
        
        # Construct base path for this sample's files
        base = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing"
        
        # Lists to track which files or steps are missing or failed
        failed_steps = []

        # Check existence and size of each expected QC file
        for qc_step in output_files_dict.keys():
            for path in output_files_dict[qc_step]:
                full_path = os.path.join(base, path.format(sid=row.Sample_ID))
                if not os.path.exists(full_path) or os.path.getsize(full_path) == 0:
                    failed_steps.append(qc_step)
        
        # If any QC step failed, log it and add sample to bad_samples
        if failed_steps:
            print(f"QC FAIL: {row.Sample_ID} from {row.flowcell_id} is missing files for QC steps: {set(failed_steps)}")
            bad_samples.add(row.Sample_ID)
    
    
    # Summary of all failures after processing
    if bad_samples:
        print(f"QC FAIL SUMMARY: {len(bad_samples)} samples missing QC files but are not in the exclude list.")
    else:
        print("QC SUCCESS! All samples have their QC files.")



##############################################################################################################################################################################

def samtools_flagstat_concat(metadata_filepath, exclude_sampleids_list, flowcells_directory, genome_shorthand, output_filepath, dtype_mapping):
    """
    Aggregates 'samtools flagstat' QC results from multiple samples into a single DataFrame and writes to a CSV.

    Parameters:
    - metadata_filepath (str) : Path to the CSV containing metadata for all samples.
    - exclude_sampleids_list (List[str]) : List of sample IDs to exclude from processing.
    - flowcells_directory (str) : Root directory where flowcell data is stored.
    - genome_shorthand (str) : Identifier for the genome reference used (e.g., 'GRCr8').
    - output_filepath (str) : Path to the CSV file where the concatenated QC results will be saved.
    - dtype_mapping (Dict[str, type]) : Dictionary specifying the data types for reading the metadata CSV.
    """
    
    # Load the multiflowcell CSV into a DataFrame
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Dictionary to hold QC metrics for all samples
    full_flagstat_dict = {}
    
    # Iterate through each row/sample in the metadata
    for _, row in metadata_df.iterrows():
        
        # Dictionary to hold metrics for a single sample
        sample_flagstat_dict = {}
        
        # Skip sample if it's in the exclusion list
        if row.Sample_ID in exclude_sampleids_list:
            continue
            
        # Construct full path to this sample's flagstat file
        filepath = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing/samtools_flagstat/{row.Sample_ID}.flagstat"
        
        # Skip sample if the file doesn't exist or is empty
        if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
            print(f"QC WARN: {row.Sample_ID} is missing flagstat file. Skipping.")
            continue
        
        # Read the flagstat file with '+' separator and no header
        df = pd.read_csv(filepath, header=None, sep="+")
        
        # Parse and compute various flagstat QC metrics
        sample_flagstat_dict["total_read_count"] =  df.iloc[0, 0]
        sample_flagstat_dict["primary_read_count"] =  df.iloc[1, 0]
        sample_flagstat_dict["primary_perc_total"] = (df.iloc[1, 0] / df.iloc[0, 0]) * 100 if pd.notnull(df.iloc[0, 0]) and df.iloc[0, 0] != 0 else float('nan')
        sample_flagstat_dict["secondary_read_count"] =  df.iloc[2, 0]
        sample_flagstat_dict["secondary_perc_total"] = (df.iloc[2, 0] / df.iloc[0, 0]) * 100 if pd.notnull(df.iloc[0, 0]) and df.iloc[0, 0] != 0 else float('nan')
        sample_flagstat_dict["supplementary_read_count"] =  df.iloc[3, 0]
        sample_flagstat_dict["duplicates_read_count"] =  df.iloc[4, 0]
        sample_flagstat_dict["primary_duplicates_read_count"] =  df.iloc[5, 0]
        sample_flagstat_dict['primary_duplicates_perc_primary'] = (df.iloc[5, 0] / df.iloc[1, 0]) * 100 if pd.notnull(df.iloc[1, 0]) and df.iloc[1, 0] != 0 else float('nan')
        sample_flagstat_dict["mapped_read_count"] =  df.iloc[6, 0]
        sample_flagstat_dict["primary_mapped_read_count"] =  df.iloc[7, 0]
        sample_flagstat_dict['primary_mapped_perc_primary'] = (df.iloc[7, 0] / df.iloc[1, 0]) * 100 if pd.notnull(df.iloc[1, 0]) and df.iloc[1, 0] != 0 else float('nan')
        sample_flagstat_dict["paired_in_sequencing_count"] =  df.iloc[8, 0]
        sample_flagstat_dict["read1_count"] =  df.iloc[9, 0]
        sample_flagstat_dict["read2_count"] =  df.iloc[10, 0]
        sample_flagstat_dict["properly_paired_read_count"] =  df.iloc[11, 0]
        sample_flagstat_dict["with_itself_and_mate_mapped_read_count"] =  df.iloc[12, 0]
        sample_flagstat_dict['with_itself_and_mate_mapped_perc_primary'] = (df.iloc[12, 0] / df.iloc[1, 0]) * 100 if pd.notnull(df.iloc[1, 0]) and df.iloc[1, 0] != 0 else float('nan')
        sample_flagstat_dict["singletons_count"] =  df.iloc[13, 0]
        sample_flagstat_dict["with_mate_mapped_diff_chr_read_count"] =  df.iloc[14, 0]
        sample_flagstat_dict["with_mate_mapped_diff_chr_mapq5_read_count"] =  df.iloc[15, 0]
        sample_flagstat_dict["with_mate_mapped_diff_chr_mapq5_perc_primary"] =  (df.iloc[15, 0] / df.iloc[1, 0]) * 100 if pd.notnull(df.iloc[1, 0]) and df.iloc[1, 0] != 0 else float('nan')
        
        # Store the metrics under the sample's ID
        full_flagstat_dict[row.Sample_ID] = sample_flagstat_dict
        
    # Convert the dictionary to a DataFrame, using Sample_ID as a column
    flagstat_ouput_df = pd.DataFrame(full_flagstat_dict).T.reset_index(names="Sample_ID")
    
    # Save the final dataframe to CSV
    flagstat_ouput_df.to_csv(output_filepath)
    
    # Return the final dataframe
    return flagstat_ouput_df


##############################################################################################################################################################################

def samtools_idxstats_concat(metadata_filepath, exclude_sampleids_list, flowcells_directory, genome_shorthand, output_filepath, dtype_mapping):
    """
    Aggregates 'samtools idxstats' QC results from multiple samples into a single DataFrame and writes to a CSV.

    Parameters:
    - metadata_filepath (str) : Path to the CSV containing metadata for all samples.
    - exclude_sampleids_list (List[str]) : List of sample IDs to exclude from processing.
    - flowcells_directory (str) : Root directory where flowcell data is stored.
    - genome_shorthand (str) : Identifier for the genome reference used (e.g., 'GRCr8').
    - output_filepath (str) : Path to the CSV file where the concatenated QC results will be saved.
    - dtype_mapping (Dict[str, type]) : Dictionary specifying the data types for reading the metadata CSV.
    """

    # Load the multiflowcell CSV into a DataFrame
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # Define chromosomes of interest
    chrom_list = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
        "chrX", "chrY", "chrM"
    ]

    # List to hold reshaped idxstats dataframes for each sample
    full_idxstats_list = []
    
    # Iterate over each row/sample in the metadata
    for _, row in metadata_df.iterrows():
        
        # Skip excluded samples
        if row.Sample_ID in exclude_sampleids_list:
            continue
            
        # Construct path to the sample's idxstats file
        filepath = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing/samtools_idxstats/{row.Sample_ID}.idxstats"
        
        # Skip sample if the file doesn't exist or is empty
        if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
            print(f"QC WARN: {row.Sample_ID} is missing idxstats file. Skipping.")
            continue
        
        # Read the idxstats file with tab separator, no header, and set custom column names
        tempdf = pd.read_csv(filepath, sep="\t", header=None, 
                             names=["chrom", "sequence_length", "count_mapped_read_segments", "count_unmapped_read_segments"])

        # Keep only standard chromosomes
        tempdf = tempdf[tempdf.chrom.isin(chrom_list)]
        
        # Reshape the table into long format with (chrom, metric) as index and sample ID as the value column
        tempdf_melted = pd.melt(tempdf, id_vars='chrom', var_name='idxstat_name', value_name=f'{row.Sample_ID}').set_index(['chrom', 'idxstat_name'])

        # Append this reshaped DataFrame to the list
        full_idxstats_list.append(tempdf_melted)
    
    # Concatenate all melted sample DataFrames along the column axis
    full_idxstats_df = pd.concat(full_idxstats_list, axis=1, join='outer').reset_index()
    
    # Save the final dataframe to a CSV
    full_idxstats_df.to_csv(output_filepath)

    # Return the final dataframe
    return full_idxstats_df


##############################################################################################################################################################################

def mosdepth_concat(metadata_filepath, exclude_sampleids_list, flowcells_directory, genome_shorthand, output_directory, dtype_mapping):
    """
    Aggregates 'mosdepth' QC results from multiple samples into a single DataFrame and writes to a CSV.

    Parameters:
    - metadata_filepath (str) : Path to the CSV containing metadata for all samples.
    - exclude_sampleids_list (List[str]) : List of sample IDs to exclude from processing.
    - flowcells_directory (str) : Root directory where flowcell data is stored.
    - genome_shorthand (str) : Identifier for the genome reference used (e.g., 'GRCr8').
    - output_filepath (str) : Path to the CSV file where the concatenated QC results will be saved.
    - dtype_mapping (Dict[str, type]) : Dictionary specifying the data types for reading the metadata CSV.
    """
    
    # Load the multiflowcell CSV into a DataFrame
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Define chromosomes of interest
    chrom_list = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
        "chrX", "chrY", "chrM"
    ]

    input_filepath_suffixes = [f".mosdepth.summary.txt", f"_excluded_duplicates.mosdepth.summary.txt", f"_excluded_unmapped.mosdepth.summary.txt"]
    output_filepath_suffixes = [f"classic", f"excluded_duplicates", f"excluded_unmapped"]
    

    for input_filepath_suffix, output_filepath_suffix in zip(input_filepath_suffixes, output_filepath_suffixes):

        # List to store melted data from all samples
        full_mosdepth_list = []
        
        # Loop through each sample row in metadata
        for _, row in metadata_df.iterrows():
            
            # Skip samples in exclusion list
            if row.Sample_ID in exclude_sampleids_list:
                continue
            
            # Construct path to the sample's mosdepth summary file
            # input_filepath = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing/mosdepth/{row.Sample_ID}.mosdepth.summary.txt"
            input_filepath = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing/mosdepth/{row.Sample_ID}{input_filepath_suffix}"
            output_filepath = f"{output_directory}/mosdepth_concat_{output_filepath_suffix}.csv"            
            
            # Skip if the file is missing or empty
            if not os.path.exists(input_filepath) or os.path.getsize(input_filepath) == 0:
                print(f"QC WARN: {row.Sample_ID} is missing mosdepth file. Skipping.")
                continue

            # Load mosdepth summary into a DataFrame
            tempdf = pd.read_csv(input_filepath, sep="\t")
            
            # Keep only rows corresponding to standard chromosomes
            tempdf = tempdf[tempdf.chrom.isin(chrom_list)]

            # Reshape to long format: (chrom, stat name) as index, sample ID as value column
            tempdf_melted = pd.melt(tempdf, id_vars='chrom', var_name='mosdepth_stat_name', value_name=f'{row.Sample_ID}').set_index(['chrom', 'mosdepth_stat_name'])
        
            # Append reshaped data for this sample
            full_mosdepth_list.append(tempdf_melted)
        
        # Combine all samples’ melted data into one wide-format DataFrame
        full_mosdepth_df = pd.concat(full_mosdepth_list, axis=1, join='outer').reset_index()
        
        # Save the final dataframe to CSV
        full_mosdepth_df.to_csv(output_filepath)

        # Return the final dataframe
        # return full_mosdepth_df

##############################################################################################################################################################################

def samtools_depth_concat(metadata_filepath, exclude_sampleids_list, flowcells_directory, genome_shorthand, dtype_mapping, output_directory):
    """
    Aggregates samtools depth QC results from multiple samples and multiple chromosomes into a single DataFrame and writes to a CSV.

    Parameters:
    - metadata_filepath (str) : Path to the CSV containing metadata for all samples.
    - exclude_sampleids_list (List[str]) : List of sample IDs to exclude from processing.
    - flowcells_directory (str) : Root directory where flowcell data is stored.
    - genome_shorthand (str) : Identifier for the genome reference used (e.g., 'GRCr8').
    - output_filepath (str) : Path to the CSV file where the concatenated QC results will be saved.
    - dtype_mapping (Dict[str, type]) : Dictionary specifying the data types for reading the metadata CSV.
    """
    
    chrom_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                 "X", "Y", "M"]
    
    # Load the multiflowcell CSV into a DataFrame
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)[1:2]

    for chrom in chrom_list:
        
        variant_coverage_dict = {}
        
        # Iterate through each sample record
        for _, row in metadata_df.iterrows():

            # Skip excluded samples
            if row.Sample_ID in exclude_sampleids_list:
                continue

            # Construct path to sample's .depth file
            filepath = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing/samtools_depth/{row.Sample_ID}_chr{chrom}.depth"

            # Skip missing or empty files and warn user
            if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
                print(f"QC WARN: {row.Sample_ID} is missing samtools depth file. Skipping.")
                continue

            print(f"Processing {row.Sample_ID}: {filepath}")
            # Read the .depth file (tab-separated): chrom, pos, depth
            tempdf = pd.read_csv(filepath, sep="\t", names=['chrom', 'pos', 'coverage_value'])

            # Add this sample's data to the list
            for _, row2 in tempdf.iterrows():
                coverage_value = row2['coverage_value']
                keyval = f"{row2.chrom}_{row2.pos}"
                if keyval not in variant_coverage_dict:
                    variant_coverage_dict[keyval] = coverage_value

                else:
                    variant_coverage_dict[keyval] += coverage_value

        # Concatenate all sample DataFrames by columns (same index: chrom, pos)
        variant_coverage_list = [{'snp': keyval, 'count': count} for keyval, count in variant_coverage_dict.items()]
        variant_coverage_df = pd.DataFrame(variant_coverage_list)

        # Save the final dataframe to CSV
        variant_coverage_df.to_csv(f"{output_directory}/samtools_depth_concat_chr{chrom}.csv", index=False)


##############################################################################################################################################################################

def sex_svm_concat(metadata_filepath, exclude_sampleids_list, genome_shorthand, flowcells_directory, dtype_mapping, output_filepath):
    
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    sex_svm_dict = {}
    
    # Iterate through each sample record
    for _, row in metadata_df.iterrows():

        # Skip excluded samples
        if row.Sample_ID in exclude_sampleids_list:
            continue

        # Construct path to sample's .depth file
        templist = ["chrX", "chrY", "allchrom"]
        
        for temp in templist:
            filepath = f"{flowcells_directory}/{row.flowcell_id}/{genome_shorthand}/03_sample_processing/samtools_view/{row.Sample_ID}_{temp}_count.csv"
            if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
                # print(f"QC WARN: {row.Sample_ID} is missing samtools depth file {filepath}. Skipping.")
                continue

            with open(filepath) as f:
                read_count = f.readline().strip()

                if row.Sample_ID not in sex_svm_dict:
                    sex_svm_dict[row.Sample_ID] = {temp: read_count}
                else:
                    sex_svm_dict[row.Sample_ID].update({temp: read_count})
    
    sex_svm_df = pd.DataFrame(sex_svm_dict).T.reset_index(names="Sample_ID")
    
    # Save the final dataframe to CSV
    sex_svm_df.to_csv(output_filepath, header=True, index=False)
    
    return sex_svm_df