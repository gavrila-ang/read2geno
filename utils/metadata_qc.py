import pandas as pd
import sys
import os
import string
import shutil
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from utils.validate_paths import validate_paths
from utils.print_aesthetics import *


def essential_columns(metadata_filepath, dtype_mapping):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    essential_columns = ['rfid', 'library_name', 'project_name', 'flowcell_id',
                         'barcode', 'pcr_barcode', 'organism', 'strain', 'sex', 
                         'fastq_files', 'seq_method', 'tissue_type']

    missing = [col for col in essential_columns if col not in metadata_df.columns]

    if missing:
        print(f"{FAIL} QC FAIL: missing essential columns: {missing}.")
        qc_pass = False
    else:
        print(f"{PASS} QC PASS: all essential columns present.")
        qc_pass = True
    
    return qc_pass


def rfid_multizero(metadata_filepath, dtype_mapping):
    
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    bad_rfids = metadata_df[metadata_df.rfid=="933000000000000"].index.tolist()

    if bad_rfids:
        print(f"{FAIL} QC FAIL: {len(bad_rfids)} samples labelled as 933000000000000.")
        qc_pass = False
    else:
        print(f"{PASS} QC PASS: no samples labelled as 933000000000000.")
        qc_pass = True
            
    return qc_pass


def rfid_scinotation(metadata_filepath, dtype_mapping):
    
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    bad_rfids = metadata_df[metadata_df.rfid=="9.33E+14"].index.tolist()

    if bad_rfids:
        print(f"{FAIL} QC FAIL: {len(bad_rfids)} samples labeled as 9.33E+14.")
        qc_pass = False
    else:
        print(f"{PASS} QC PASS: no samples labeled as 9.33E+14.")
        qc_pass = True
            
    return qc_pass


def rfid_short(metadata_filepath, dtype_mapping, current_round, log_dir, len_limit=4):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    bad_rfids = metadata_df.loc[metadata_df["rfid"].apply(lambda x: len(str(x).strip()) < 8)]

    if bad_rfids.shape[0] != 0:
        print(f"{WARN} QC WARN: {len(bad_rfids)} samples with rfids less than {len_limit} in length.")
        bad_rfids[['rfid', 'Library_ID', 'Sample_Project']].to_csv(f"{log_dir}/{current_round}_short_rfids.csv", index=None)        
        print(f"{FILE} Saving list of short rfids to {log_dir}/{current_round}_short_rfids.csv")
        qc_pass = False        

    else:
        print(f"{PASS} QC PASS: no samples with rfids less than {len_limit} in length. Continuing metadata QC...")
        qc_pass = True
        
    return qc_pass


def rfid_missing_data(metadata_filepath, dtype_mapping, current_round, log_dir):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    required_fields = [
        "rfid", "Library_ID", "Sample_Project", "flowcell_id",
        "Sample_Barcode", "organism", "strain", "sex", "Fastq_Files",
        "seq_method"
    ]

    has_rfid_errors = 0

    # lists, dict to store metadata for problem samples
    lib_names = []
    proj_names = []
    rfids = []
    fields = []
    missing_fields = {}

    for _, row in metadata_df.iterrows():

        # error : check if any rfid is not associated with a required_fields
        missing_required_fields = [field for field in required_fields if pd.isnull(row[field])]
        
        if missing_required_fields:
            # print(f"{FAIL} QC FAIL: RFID '{row.rfid}' from library '{row.library_name}' and project_name '{row.project_name}' is missing the following required fields: {missing_required_fields}")
            has_rfid_errors += 1
            lib_names.append(row['Library_ID'])
            proj_names.append(row['Sample_Project'])
            rfids.append(row.rfid)
            fields.append('|'.join(missing_required_fields))
            
            for field in missing_required_fields:
                if field not in missing_fields:
                    missing_fields[field] = []
                missing_fields[field].append(row.rfid)
        
        # warn : check if any rfid is not associated with non-required fields
        missing_optional_fields = [field for field in metadata_df.columns if (field not in required_fields) & (pd.isnull(row[field]))]
        if missing_optional_fields:
            continue        

    if has_rfid_errors != 0:
        qc_pass = False
        
        # log a csv of samples with missing data
        df = pd.DataFrame({
            'Library_ID': lib_names,
            'rfid': rfids,
            'Sample_Project': proj_names,
            'missing_fields': fields
        })
        outfile = os.path.join(log_dir, f'{current_round}_qc_fail_rfids_missing_data.csv')
        df.to_csv(outfile, index=False)

        # save files of rfids for each column with missing data
        for field, rfids in missing_fields.items():
            df = pd.DataFrame({field: rfids})
            field_outfile = os.path.join(log_dir, f"{current_round}_qc_fail_rfids_missing_{field}.csv")
            df.to_csv(field_outfile, index=False, header=False)

        # print(f"{FAIL} QC FAIL: {has_rfid_errors} samples with rfids missing essential data.")

        # count missing rfids per missing field
        for field, rfids in missing_fields.items():
            count = len(rfids)
            if field != 'sex':
                print(f"{FAIL} QC FAIL: {count} samples are missing {field} data.")
            else:
                print(f"{WARN} QC WARN: {count} samples are missing {field} data.")

        print(f"{FILE} {has_rfid_errors} samples with missing essential data written to {outfile}")

    else:
        print(f"{PASS} QC PASS: all samples have essential data. Continuing metadata QC...")
        qc_pass = True
        
    return qc_pass


def fastq_validity(metadata_filepath, dtype_mapping, sequencing_directory, current_round, log_dir):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    metadata_df['malformed_fastq'] = 0
    metadata_df['missing_fastq'] = 0
    metadata_df['fastq_1'] = None
    metadata_df['fastq_2'] = None
    malformed_fastqs = 0
    missing_fastqs = 0

    for i, row in metadata_df.iterrows():
        
        # error : check if any rfid has a malformed fastq_file
        if row.seq_method in ['riptide', 'LSW_lcwgs']:
            try:
                fq1_filename, fq2_filename = row["Fastq_Files"].split(";")
                fq1_filepath = os.path.join(sequencing_directory, row["flowcell_id"], fq1_filename.strip())
                fq2_filepath = os.path.join(sequencing_directory, row["flowcell_id"], fq2_filename.strip())

                if not os.path.exists(fq1_filepath):
                    missing_fastqs += 1
                    metadata_df.loc[i, 'missing_fastq'] = 1
                    metadata_df.loc[i, 'fastq_1'] = fq1_filepath

                if not os.path.exists(fq2_filepath):
                    missing_fastqs += 1
                    metadata_df.loc[i, 'missing_fastq'] = 1
                    metadata_df.loc[i, 'fastq_2'] = fq2_filepath

            except ValueError:
                malformed_fastqs += 1
                metadata_df.loc[i, 'malformed_fastq'] = 1
                print(f"{FAIL} QC FAIL: RFID '{row.rfid}' from library '{row.library_name}' and seq_method '{row.seq_method}' has a malformed fastq_files field: '{row['fastq_files']}'")
                
        elif row.seq_method == 'gbs':
            try:
                fq1_filepath = os.path.join(sequencing_directory, row["flowcell_id"], row["Fastq_Files"].strip())
                
                if not os.path.exists(fq1_filepath):
                    missing_fastqs += 1
                    metadata_df.loc[i, 'missing_fastq'] = 1
                    metadata_df.loc[i, 'fastq_1'] = fq1_filepath
                    metadata_df.loc[i, 'fastq_2'] = "NA"                    

            except ValueError:
                malformed_fastqs += 1
                metadata_df.loc[i, 'malformed_fastq'] = 1     

        else:
            print(f"Invalid sequencing method for {row.Sample_ID}: {row.seq_method}")         

    if malformed_fastqs != 0:
        qc_pass = False
        malformed_df = metadata_df[metadata_df['malformed_fastq']==1]
        malformed_df = malformed_df[['rfid','Library_ID','Sample_Project','fastq_1','fastq_2']]
        outfile = os.path.join(log_dir, f'{current_round}_qc_fail_fastq_malformed_files.csv')
        malformed_df.to_csv(outfile, index=False)
        print(f'{FAIL} QC FAIL: {malformed_df.shape[0]} samples with malformed fastq file paths written to {outfile}')
    
    if missing_fastqs != 0:
        qc_pass = False
        missing_df = metadata_df[metadata_df['missing_fastq']==1]
        missing_df = missing_df[['rfid','Library_ID','Sample_Project','fastq_1','fastq_2']]
        outfile = os.path.join(log_dir, f'{current_round}_qc_fail_fastq_missing_files.csv')
        missing_df.to_csv(outfile, index=False)
        print(f'{FAIL} QC FAIL: {missing_df.shape[0]} samples with missing fastq file paths written to {outfile}')

    if (malformed_fastqs == 0 and missing_fastqs == 0):
        print(f"{PASS} QC PASS: ALL samples have valid fastq filepaths. Continuing metadata QC...")
        qc_pass = True
        
    return qc_pass


def library_association(metadata_filepath, dtype_mapping):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    bad_libraries = 0
    
    for library, group in metadata_df.groupby("Library_ID"):
        
        flowcell_ids = group["flowcell_id"].unique()
        if len(flowcell_ids) > 1:
            print(f"{WARN} WARNING: Library '{library}' is associated with multiple flowcell_ids: {flowcell_ids}")

        fastq_files = group["Fastq_Files"].unique()
        if len(fastq_files) > 1:
            print(f"{WARN} WARNING: Library '{library}' is associated with multiple fastq_files: {fastq_files}")

        seq_methods = group["seq_method"].unique()
        if len(seq_methods) > 1:
            print(f"{WARN} WARNING: Library '{library}' is associated with multiple seq_method: {seq_methods}")

        if (len(flowcell_ids) > 1) or (len(fastq_files) > 1) or (len(seq_methods) > 1):
            bad_libraries += 1

    if bad_libraries != 0:
        print(f"{FAIL} QC FAIL: {bad_libraries} libraries associated with too many variables. Exiting pipeline.")
        qc_pass = False
    else:
        print(f"{PASS} QC PASS: each library is associated with only one flowcell_id, one set of fastq_files, seq_method. Continuing metadata QC...")
        qc_pass = True


def fastq_file_association(metadata_filepath, dtype_mapping):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    bad_fastq_files = 0
    
    for fastq_files, group in metadata_df.groupby("Fastq_Files"):

        libraries = group["Library_ID"].unique()

        if len(libraries) > 1:
            print(f"{WARN} WARNING: fastq_files set '{fastq_files}' is associated with multiple libraries: {libraries}")
            bad_fastq_files += 1

    if bad_fastq_files != 0:
        print(f"{FAIL} QC FAIL: {bad_fastq_files} fastq_files associated with too many libraries. Exiting pipeline.")
        qc_pass = False
        
    else:
        print(f"{PASS} QC PASS: each set of fastq_files is associated with only one library. Continuing metadata QC...")
        qc_pass = True


def rename_columns(metadata_filepath, dtype_mapping, output_filepath):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # Rename columns
    # "rfid": "Sample_Name"
    rename_col_dict = {
        "library_name": "Library_ID",  
        "project_name": "Sample_Project",
        "barcode": "Sample_Barcode",
        "fastq_files": "Fastq_Files",
        "pcr_barcode": "PCR_Barcode"
    }
    
    metadata_df = metadata_df.rename(columns=rename_col_dict)

    # Save renamed columns 
    metadata_df.to_csv(output_filepath, index=False)
    print(f"{FILE} Renamed metadata columns. Saving all new flowcell sheets to {output_filepath}.")


def create_sampleid(metadata_filepath, dtype_mapping, output_filepath):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # Create Sample_ID column

    metadata_df["Sample_ID"] = metadata_df.apply(lambda row: f"{row['Library_ID']}_{row['rfid']}_{row['Sample_Barcode']}" \
                                                                            if row['seq_method'] in ['riptide', 'LSW_lcwgs'] \
                                                                            else f"{row['Library_ID']}_{row['rfid']}", axis=1)

    # Create Sample_Name column
    metadata_df['Sample_Name'] = metadata_df.apply(lambda row: f"{row['rfid']}_{row['Sample_Barcode']}" \
                                                                            if row['seq_method'] in ['riptide', 'LSW_lcwgs'] \
                                                                            else f"{row['rfid']}", axis=1)

    # Reorder columns with Sample_ID at front
    identifier_cols = ["Sample_ID", "Sample_Name", "rfid", "Sample_Barcode", "Library_ID", "flowcell_id"]
    new_col_order = identifier_cols + [col for col in metadata_df.columns if col not in identifier_cols]
    metadata_df = metadata_df[new_col_order]

    # Save renamed metadata 
    metadata_df.to_csv(output_filepath, index=False)
    print(f"{FILE} Created Sample_ID column. Saving all new flowcell sheets to {output_filepath}.")


def sampleid_missing_data(metadata_filepath, exclude_sampleids_list, dtype_mapping):
    
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    required_fields = [
        "rfid", "Library_ID", "Sample_Project", "flowcell_id",
        "organism", "strain", "sex",
        "seq_method"
    ]

    has_sampleid_errors = 0
    
    for _, row in metadata_df.iterrows():
        
        if row.Sample_ID not in exclude_sampleids_list:

            # error : check if any rfid is not associated with a required_fields
            missing_required_fields = [field for field in required_fields if pd.isnull(row[field])]
            if missing_required_fields:
                print(f"{FAIL} QC FAIL: Sample ID '{row.Sample_ID}' from library '{row.Library_ID}' and project_name '{row.Sample_Project}' is missing the following required fields: {missing_required_fields}")
                has_sampleid_errors += 1
            # warn : check if any rfid is not associated with non-required fields
            missing_optional_fields = [field for field in metadata_df.columns if (field not in required_fields) & (pd.isnull(row[field]))]
            if missing_optional_fields:
                continue        

    if has_sampleid_errors != 0:
        print(f"{FAIL} QC FAIL: {has_sampleid_errors} samples with rfids missing essential data. Exiting pipeline.")
        qc_pass = False
        
    else:
        print(f"{PASS} QC PASS: all samples have essential data. Continuing metadata QC...")
        qc_pass = True
        
    return qc_pass


def hsrats_only(metadata_filepath, dtype_mapping, output_filepath):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Filter for hs rat samples
    metadata_df["strain"] = metadata_df["strain"].str.lower()
    metadata_df["organism"] = metadata_df["organism"].str.lower()
    filtered = metadata_df[(metadata_df.strain=="heterogeneous stock") & (metadata_df.organism=="rat")]

    # Save filtered allnew_fs 
    filtered.to_csv(output_filepath, index=False)
    print(f"{FILE} Filtered for heterogenous stock rats only. Saving all new flowcell sheets to {output_filepath}.")
    
    # Notify user if there are non-hsrat samples
    if len(filtered) < len(metadata_df):
        count_non_hsrat = len(metadata_df) - len(filtered)
        print(f"{WARN} WARNING: There are {count_non_hsrat} non-HS rat samples in {metadata_filepath}.")
    else:
        print(f"{PASS} QC PASS: All samples are HS rat samples in {metadata_filepath}.")


def remove_bad_sampleids(metadata_filepath, current_round, drop_sample_log_filepath, drop_library_log_filepath, drop_project_log_filepath, reinclude_sample_log_filepath, dtype_mapping, output_filepath):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    dropped_samples = pd.read_csv(drop_sample_log_filepath).query(f'round_dropped != @current_round').Sample_ID.tolist()
    dropped_libraries = pd.read_csv(drop_library_log_filepath).query(f'round_dropped != @current_round').library_name.tolist()
    dropped_projects = pd.read_csv(drop_project_log_filepath).query(f'round_dropped != @current_round').project_name.tolist()    
    reinclude_samples = pd.read_csv(reinclude_sample_log_filepath).Sample_ID.tolist()

    current_drops_log = metadata_df.query('Sample_ID.isin(@dropped_samples) or \
                                    Library_ID.isin(@dropped_libraries) or \
                                    Sample_Project.isin(@dropped_projects)'
                                   )
    
    # Drop any sample that's in the reinclude log
    current_drops_log = current_drops_log.query('~Sample_ID.isin(@reinclude_samples)')
    current_exclude_sampleids_list = current_drops_log.Sample_ID.tolist()
    filtered_df = metadata_df[~metadata_df.Sample_ID.isin(current_exclude_sampleids_list)]    

    filtered_df.to_csv(output_filepath, index=False)
    print(f"{FILE} Saving filtered multi-flowcell sheet to {output_filepath}.")

    if len(current_exclude_sampleids_list) != 0:
        print(f"There are {len(current_exclude_sampleids_list)} samples currently in the exclude sample IDs list.")
        print(f"There are {filtered_df.shape[0]} remaining samples in the multi-flowcell sheet after removal of unwanted Sample IDs.")


def detect_dup_sampleids(metadata_filepath, exclude_sampleids_list, dtype_mapping):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    filtered_df = metadata_df[~metadata_df.Sample_ID.isin(exclude_sampleids_list)]

    # Identify duplicated Sample_IDs
    duplicate_sample_ids = filtered_df.loc[filtered_df.Sample_ID.duplicated(keep=False),
                                            ["Sample_ID", "Sample_Name", "Library_ID", "flowcell_id"]
                                            ].copy().sort_values(by="Sample_Name")

    if len(duplicate_sample_ids != 0):
        print(f"{FAIL} QC FAIL: There are {len(duplicate_sample_ids)} duplicate Sample_IDs in the dataset.")
        print(duplicate_sample_ids)
        qc_pass = False

    else:
        print(f"{PASS} QC PASS: There are no duplicate Sample_IDs in the dataset.")
        qc_pass = True

    return qc_pass


def rename_dup_rfids(metadata_filepath, dtype_mapping, output_filepath, output_summary_filepath):

    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # Identify duplicated Sample_IDs
    duplicate_sample_ids = metadata_df.loc[metadata_df.Sample_ID.duplicated(keep=False),
                                                    ["Sample_Name", "Sample_ID", "Sample_Barcode", "tissue_type", "Library_ID", "flowcell_id", "Sample_Project"]
                                                    ].copy().sort_values(by="Sample_Name")

    # Generate alpha suffixes for duplicates; assigns a suffix like _26, _27, etc., when x > 26
    duplicate_sample_ids["suffix"] = duplicate_sample_ids.groupby("Sample_ID").cumcount().map(
        lambda x: string.ascii_lowercase[x] if x < 26 else f'_{x}')

    # Generates per-sample duplicate counts
    sampleid_duplicate_counts = duplicate_sample_ids['Sample_ID'].value_counts(dropna=False)
    duplicate_sample_ids = duplicate_sample_ids.merge(sampleid_duplicate_counts, on='Sample_ID', how='left')

    # Create modified Sample_IDs and Sample_Names
    duplicate_sample_ids["modified_Sample_ID"] = duplicate_sample_ids["Sample_ID"] + duplicate_sample_ids["suffix"]
    duplicate_sample_ids["modified_Sample_Name"] = duplicate_sample_ids["Sample_Name"] + duplicate_sample_ids["suffix"]

    # Retain original names
    metadata_df['Sample_ID_original'] = metadata_df['Sample_ID'].copy()
    metadata_df['Sample_Name_original'] = metadata_df['Sample_Name'].copy()

    # Perform replacements
    for _, row in duplicate_sample_ids.iterrows():
        condition = (
            (metadata_df['Sample_ID'] == row['Sample_ID']) &
            (metadata_df['Sample_Barcode'] == row['Sample_Barcode']) &
            (metadata_df['tissue_type'] == row['tissue_type']) &
            (metadata_df['flowcell_id'] == row['flowcell_id']) &
            (metadata_df['Sample_Project'] == row['Sample_Project'])
        )
        try:
            metadata_df.loc[condition, 'Sample_ID'] = row['modified_Sample_ID']
            metadata_df.loc[condition, 'Sample_Name'] = row['modified_Sample_Name']
        except Exception as e:
            print(e)

    # Reorder columns
    new_col_order = ["Sample_ID", "Sample_ID_original", "Sample_Name", "Sample_Name_original"] + \
                    [col for col in metadata_df.columns if col not in ["Sample_ID", "Sample_ID_original", "Sample_Name", "Sample_Name_original"]]
    metadata_df = metadata_df[new_col_order]

    # Print how many duplicated sample_ids
    if len(duplicate_sample_ids) != 0:
        print(f"{WARN} WARNING: There are {len(duplicate_sample_ids)} Sample_IDs that are duplicated. Will append alpha characters to Sample_ID and Sample_Name for fgbio-compatibility.")

    # Verify no duplicated Sample_IDs remain
    remaining_dups = metadata_df["Sample_ID"].duplicated().sum()
    if remaining_dups > 0:
        print(f"{FAIL} QC FAIL: {remaining_dups} duplicate Sample_IDs still remain in the dataset. Check output at {output_filepath}. Exiting.")
        qc_pass = False
    else:
        print(f"{PASS} QC PASS: No duplicate Sample_IDs remaining in the dataset.")
        qc_pass = True

    # Save updated dataframe and summary of changes
    duplicate_sample_ids.to_csv(output_summary_filepath, index=False)
    metadata_df.to_csv(output_filepath, index=False)
    print(f"{FILE} Saving all new flowcell sheets to {output_filepath}.")
    print(f"{FILE} Saving summary to {output_summary_filepath}.")

    return qc_pass


def make_library_sheets(metadata_filepath, exclude_sampleids_list, dtype_mapping, output_dir):
    
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    filtered_df = metadata_df[~metadata_df.Sample_ID.isin(exclude_sampleids_list)]

    total_rows = 0
    for library_flowcell, group in filtered_df.groupby(["Library_ID", "flowcell_id"]):

        library_id = library_flowcell[0]
        flowcell_id = library_flowcell[1]

        lib_sheet = group[["Sample_ID", "Sample_Name", "Library_ID", "Sample_Barcode"]]
        total_rows += len(lib_sheet)

        output_filepath = os.path.join(output_dir, f"{flowcell_id}_{library_id}_library_sheet.csv")
        lib_sheet.to_csv(output_filepath, index=False)
        # print(f"Created library-level sheet for {library_flowcell[0]}. Saved to {output_filepath}.")

    # Final check
    unassigned_sampleids = len(filtered_df) - total_rows
    if unassigned_sampleids != 0:
        print(f"{FAIL} QC ERROR: {unassigned_sampleids} Sample_IDs have not been allocated to a library-specific sheet.")
        qc_pass = False
    else:
        print(f"{PASS} QC PASS: All Sample_IDs have been allocated to a library-specific sheet.")
        qc_pass = True
        
    return qc_pass



