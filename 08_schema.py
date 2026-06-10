#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 8: Generate database-ready QC log files
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Cleaned metadata file
    - Drop library/project logs
    - Alignment files (samtools flagstat, mosdepth)
    - Sex QC SVM output files
    - PLINK heterozygosity and missingness files (.het, .smiss)

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Drop library log (.csv)
    - Drop project log (.csv)
    - Metadata QC log (.csv)
    - Demultiplexing QC log (.csv)
    - BAM QC log (.csv)
    - Sex QC log (.csv)
    - Heterozygosity QC log (.csv)
    - Missingness QC log (.csv)

Assumptions:
    - Custom input file directory structure
    - Samples have been processed through steps 1-7

"""

##############################################################################################################################################################################

# === Load python packages ===
import user
import utils
import os
import sys
import subprocess
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from tqdm import tqdm
from time import time

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Load system arguments ===
output_directory = str(sys.argv[1])
genome_size = int(sys.argv[2])

print(f"output_directory: {output_directory}")
print(f"genome_size: {genome_size}")

# === Declare input/output filepaths ===
clean_metadata_filepath = f"{round_directory}/output/01_metqc/flowcell_sheets/{current_round}_metadata_cleaned.csv"
output_prefix = f"{output_directory}/{current_round}"

# === Declare data type of special columns ===
dtype_mapping = {'rfid': str}

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
    elapsed = time() - start
    print(f"Time taken: {elapsed}")

def sh(cmd):
    """
    cmd: list[str]
    """
    subprocess.run(cmd, check=True)

##############################################################################################################################################################################

# === Functions ===

def create_drop_library_log(code, current_round, output_filepath):
    """
    Purpose: Filter and export the drop library log for the current round.
    """
    df = pd.read_csv(f"{code}/user/drop_library_log.csv")
    df = df.query('round_dropped == @current_round')
    df.to_csv(output_filepath, index=False)
    return df


def create_drop_project_log(code, current_round, output_filepath):
    """
    Purpose: Filter and export the drop project log for the current round.
    """
    df = pd.read_csv(f"{code}/user/drop_project_log.csv")
    df = df.query('round_dropped == @current_round')
    df.to_csv(output_filepath, index=False)
    return df


def create_metadata_qc(sequencing_data, multi_flowcell_filepath, current_round, output_filepath):
    """
    Purpose: Creates a database-ready table that summarises the quality checks performed on the metadata.
    The `unresolvable_sid` column remains empty, as it is a column that is manually filled in by the user
    if a sample fails a quality control check that is not listed below.
    """
    tempdf = pd.read_csv(multi_flowcell_filepath, dtype={'rfid': 'str'})

    fq_dict = {}

    for index, row in tqdm(tempdf.iterrows(), total=len(tempdf)):
        sid = row['Sample_ID']
        if row['seq_method'] == 'gbs':
            fq1_filepath = '/tscc' + row['Fastq_Files']
            library_fastq_exists = os.path.exists(fq1_filepath)
            fq_dict[sid] = {'fq1_filepath': fq1_filepath, 'fq2_filepath': None, 'library_fastq_exists': library_fastq_exists}
        else:
            fq1, fq2 = row['Fastq_Files'].split(';')
            fq1_filepath = sequencing_data + '/' + row['flowcell_id'] + '/' + fq1
            fq2_filepath = sequencing_data + '/' + row['flowcell_id'] + '/' + fq2.strip()
            library_fastq_exists = os.path.exists(fq1_filepath) and os.path.exists(fq2_filepath)
            fq_dict[sid] = {'fq1_filepath': fq1_filepath, 'fq2_filepath': fq2_filepath, 'library_fastq_exists': library_fastq_exists}

    fq_df = pd.DataFrame.from_dict(fq_dict).T
    metadata_qc_df = pd.concat([tempdf.set_index('Sample_ID'), fq_df], axis=1)

    required_fields = ["rfid", "Library_ID", "Sample_Project", "flowcell_id", "organism", "strain", "seq_method"]
    metadata_qc_df['rfid_not_missing_data_list'] = metadata_qc_df[required_fields].notna().any(axis=1)

    metadata_qc_df['strain'] = metadata_qc_df['strain'].str.lower()
    metadata_qc_df['hsrats_only'] = metadata_qc_df.apply(lambda row: row['strain'] == 'heterogeneous stock', axis=1)
    metadata_qc_df['rfid_non_multizero'] = metadata_qc_df.apply(lambda row: row['rfid'] != '933000000000000', axis=1)
    metadata_qc_df['rfid_non_scinotation'] = metadata_qc_df.apply(lambda row: row['rfid'] != '9.33E+14', axis=1)
    metadata_qc_df['rfid_not_short'] = metadata_qc_df.apply(lambda row: len(row['rfid']) > 4, axis=1)
    metadata_qc_df['metadata_qc_status'] = metadata_qc_df.apply(
        lambda row: row['rfid_not_missing_data_list'] and
                    row['hsrats_only'] and
                    row['rfid_non_multizero'] and
                    row['rfid_non_scinotation'],
        axis=1
    )

    metadata_qc_df['round_calculated'] = current_round
    metadata_qc_df['unresolvable_sid'] = False
    metadata_qc_df['notes'] = None

    metadata_qc_df = metadata_qc_df.rename(columns={'Library_ID': 'library_name', 'Sample_Barcode': 'barcode'})
    metadata_qc_df = metadata_qc_df[[
        'library_name', 'rfid', 'barcode', 'round_calculated',
        'library_fastq_exists', 'fq1_filepath', 'fq2_filepath',
        'rfid_not_missing_data_list', 'hsrats_only',
        'rfid_non_multizero', 'rfid_non_scinotation', 'rfid_not_short',
        'unresolvable_sid', 'metadata_qc_status',
        'notes'
    ]]

    metadata_qc_df = metadata_qc_df.reset_index(drop=True)
    metadata_qc_df.to_csv(output_filepath, index=False)
    print(f'Results written out to: {output_filepath}')

    return metadata_qc_df


def create_demux_qc(multi_flowcell_filepath, flowcells_directory, genome_shorthand, current_round, output_filepath):
    """
    Purpose: Creates a database-ready table that summarises the demultiplexing QC results.
    """
    tempdf = pd.read_csv(multi_flowcell_filepath, dtype={'rfid': 'str'})
    tempdf = tempdf[['Sample_ID', 'Library_ID', 'rfid', 'Sample_Barcode', 'flowcell_id', 'seq_method']].set_index('Sample_ID')

    fq_dict = {}

    for idx_val, row in tqdm(tempdf.iterrows(), total=len(tempdf)):
        sid = idx_val
        dirpath = f"{flowcells_directory}/{row['flowcell_id']}/{genome_shorthand}/02_demux/demux/{row['Library_ID']}"

        if row['seq_method'] == 'gbs':
            fq1_filepath_value = f'{dirpath}/{row["Library_ID"]}_{row["rfid"]}.fq.gz'
            fastq_exists_value = os.path.exists(fq1_filepath_value)
            fq_dict[sid] = {'fq1_filepath': fq1_filepath_value, 'fq2_filepath': 'single_read', 'fastq_exists': fastq_exists_value}
        else:
            file_prefix = f'{row["Library_ID"]}_{row["rfid"]}_{row["Sample_Barcode"]}-{row["rfid"]}_{row["Sample_Barcode"]}-{row["Sample_Barcode"]}'
            fq1_filepath_value = f"{dirpath}/{file_prefix}_R1.fastq.gz"
            fq2_filepath_value = f"{dirpath}/{file_prefix}_R2.fastq.gz"
            fastq_exists_value = os.path.exists(fq1_filepath_value) and os.path.exists(fq2_filepath_value)
            fq_dict[sid] = {'fq1_filepath': fq1_filepath_value, 'fq2_filepath': fq2_filepath_value, 'fastq_exists': fastq_exists_value}

    fq_df = pd.DataFrame.from_dict(fq_dict).T
    demux_qc_df = pd.concat([tempdf, fq_df], axis=1)

    demux_qc_df['fastq_size'] = demux_qc_df.apply(
        lambda row: os.path.getsize(row['fq1_filepath']) if row['fastq_exists'] else 0, axis=1
    )
    demux_qc_df['small_fastq_size'] = demux_qc_df.apply(lambda row: row['fastq_size'] < 10_000_000, axis=1)
    demux_qc_df['demux_qc_status'] = demux_qc_df.apply(lambda row: row['fastq_exists'], axis=1)

    demux_qc_df['round_calculated'] = current_round
    demux_qc_df['notes'] = None

    demux_qc_df = demux_qc_df.rename(columns={'Library_ID': 'library_name', 'Sample_Barcode': 'barcode'})
    demux_qc_df = demux_qc_df[[
        'library_name', 'rfid', 'barcode', 'round_calculated',
        'fq1_filepath', 'fq2_filepath',
        'fastq_exists', 'fastq_size', 'small_fastq_size',
        'demux_qc_status',
        'notes'
    ]]

    demux_qc_df = demux_qc_df.reset_index(drop=True)
    demux_qc_df.to_csv(output_filepath, index=False)
    print(f'Results written out to: {output_filepath}')

    return demux_qc_df


def create_bam_qc(multi_flowcell_filepath, flowcells_directory, genome_shorthand, genome_size, current_round, output_filepath):
    """
    Purpose: Creates a database-ready table that summarises the alignment (BAM) QC results.
    """
    tempdf = pd.read_csv(multi_flowcell_filepath, dtype={'rfid': 'str'})
    tempdf = tempdf[['Sample_ID', 'Library_ID', 'rfid', 'Sample_Barcode', 'flowcell_id']].set_index('Sample_ID')

    flagstat_dict = {}
    mosdepth_mapped_dict = {}
    mosdepth_nondup_dict = {}

    for sid, row in tqdm(tempdf.iterrows(), total=len(tempdf)):
        outdir = f'{flowcells_directory}/{row["flowcell_id"]}/{genome_shorthand}/03_sample_processing'
        samtools_flagstat = f'{outdir}/samtools_flagstat/{sid}.flagstat'
        mosdepth_mapped = f'{outdir}/mosdepth/{sid}_excluded_unmapped.mosdepth.summary.txt'
        mosdepth_nondup = f'{outdir}/mosdepth/{sid}_excluded_duplicates.mosdepth.summary.txt'

        if os.path.exists(samtools_flagstat):
            flagdf = pd.read_csv(samtools_flagstat, sep=r'[+()]', header=None, engine='python')
            flagdf.iloc[15, 1] = 'with mate mapped to a different chr mapq5'
            index_list = flagdf[1].str.findall(r"\w+").apply(lambda x: "_".join(x[1:])).tolist()
            flagdf.index = index_list
            flagstat_dict[sid] = flagdf[0]

        if os.path.exists(mosdepth_mapped):
            mosdf = pd.read_csv(mosdepth_mapped, sep='\t', index_col='chrom')
            mosdepth_mapped_dict[sid] = mosdf['bases']

        if os.path.exists(mosdepth_nondup):
            mosdf = pd.read_csv(mosdepth_nondup, sep='\t', index_col='chrom')
            mosdepth_nondup_dict[sid] = mosdf['bases']

    flagstat_concat = pd.concat(flagstat_dict).unstack()
    flagstat_concat.columns = [col + '_read_count' for col in flagstat_concat.columns]
    mosdepth_mapped_concat = pd.concat(mosdepth_mapped_dict).unstack().rename(
        columns={'total': 'total_mapped_base_count'})['total_mapped_base_count']
    mosdepth_nondup_concat = pd.concat(mosdepth_nondup_dict).unstack().rename(
        columns={'total': 'total_mapped_nonduplicate_base_count'})['total_mapped_nonduplicate_base_count']

    bam_qc_df = pd.concat([tempdf, flagstat_concat, mosdepth_mapped_concat, mosdepth_nondup_concat], axis=1)

    bam_qc_df['markdup_bam_filepath'] = (
        flowcells_directory + '/' + bam_qc_df['flowcell_id'] + '/' + genome_shorthand +
        '/03_sample_processing/picard_markdup/' + bam_qc_df.index + '_sorted_mkDup.bam'
    )
    bam_qc_df['markdup_bam_file_exists'] = bam_qc_df.apply(
        lambda row: os.path.exists(row['markdup_bam_filepath']), axis=1
    )

    bam_qc_df['round_calculated'] = current_round
    bam_qc_df['notes'] = None

    bam_qc_df['coverage'] = bam_qc_df['total_mapped_base_count'] / genome_size
    bam_qc_df['duplicate_adjusted_coverage'] = bam_qc_df['total_mapped_nonduplicate_base_count'] / genome_size
    bam_qc_df['primary_perc_total'] = (bam_qc_df['primary_read_count'] / bam_qc_df['in_total_read_count']) * 100
    bam_qc_df['secondary_perc_total'] = (bam_qc_df['secondary_read_count'] / bam_qc_df['in_total_read_count']) * 100
    bam_qc_df['primary_duplicates_perc_primary'] = (bam_qc_df['primary_duplicates_read_count'] / bam_qc_df['primary_read_count']) * 100
    bam_qc_df['primary_mapped_perc_primary'] = (bam_qc_df['primary_mapped_read_count'] / bam_qc_df['primary_read_count']) * 100
    bam_qc_df['with_itself_and_mate_mapped_perc_primary'] = (bam_qc_df['with_itself_and_mate_mapped_read_count'] / bam_qc_df['primary_read_count']) * 100
    bam_qc_df['with_mate_mapped_diff_chr_perc_primary'] = (bam_qc_df['with_mate_mapped_to_a_different_chr_read_count'] / bam_qc_df['primary_read_count']) * 100
    bam_qc_df['with_mate_mapped_diff_chr_mapq5_perc_primary'] = (bam_qc_df['mate_mapped_to_a_different_chr_mapq5_read_count'] / bam_qc_df['primary_read_count']) * 100
    bam_qc_df['primary_nonduplicates_read_count'] = bam_qc_df['primary_read_count'] - bam_qc_df['primary_duplicates_read_count']
    bam_qc_df['enough_non_duplicate_reads'] = bam_qc_df['primary_nonduplicates_read_count'] > 1

    bam_qc_df = bam_qc_df.rename(columns={
        'Library_ID': 'library_name',
        'Sample_Barcode': 'barcode',
        'in_total_read_count': 'total_read_count',
        'paired_in_sequencing_read_count': 'paired_in_sequencing_count',
        'read1_read_count': 'read1_count',
        'read2_read_count': 'read2_count',
        'singletons_read_count': 'singletons_count',
        'with_mate_mapped_to_a_different_chr_read_count': 'with_mate_mapped_diff_chr_read_count',
        'mate_mapped_to_a_different_chr_mapq5_read_count': 'with_mate_mapped_diff_chr_mapq5_read_count'
    })

    bam_qc_df['bam_qc_status'] = bam_qc_df.apply(
        lambda row: row['markdup_bam_file_exists'] and row['enough_non_duplicate_reads'], axis=1
    )

    bam_qc_df = bam_qc_df[[
        'library_name', 'rfid', 'barcode', 'round_calculated',
        'markdup_bam_filepath', 'markdup_bam_file_exists',
        'total_read_count',
        'primary_read_count', 'primary_perc_total',
        'secondary_read_count', 'secondary_perc_total',
        'supplementary_read_count',
        'duplicates_read_count', 'primary_duplicates_read_count', 'primary_duplicates_perc_primary',
        'mapped_read_count',
        'primary_mapped_read_count', 'primary_mapped_perc_primary',
        'paired_in_sequencing_count',
        'read1_count', 'read2_count',
        'properly_paired_read_count',
        'with_itself_and_mate_mapped_read_count', 'with_itself_and_mate_mapped_perc_primary',
        'singletons_count',
        'with_mate_mapped_diff_chr_read_count', 'with_mate_mapped_diff_chr_perc_primary',
        'with_mate_mapped_diff_chr_mapq5_read_count', 'with_mate_mapped_diff_chr_mapq5_perc_primary',
        'total_mapped_base_count', 'total_mapped_nonduplicate_base_count',
        'coverage', 'duplicate_adjusted_coverage',
        'primary_nonduplicates_read_count', 'enough_non_duplicate_reads',
        'bam_qc_status',
        'notes'
    ]]

    bam_qc_df = bam_qc_df.reset_index(drop=True)
    bam_qc_df.to_csv(output_filepath, index=False)
    print(f'Results written out to: {output_filepath}')

    return bam_qc_df


def create_sex_qc_log(multi_flowcell_filepath, round_directory, riptide_only, current_round,
                      riptide_for_training_threshold, riptide_undetermined_threshold,
                      gbs_for_training_threshold, gbs_undetermined_threshold, output_filepath):
    """
    Purpose: Creates a database-ready table that summarises the sex QC SVM results.
    """
    tempdf = pd.read_csv(multi_flowcell_filepath)[['Sample_ID', 'Library_ID', 'rfid', 'Sample_Barcode', 'seq_method']]

    base_path = f"{round_directory}/output/04_align_qc/sex_qc"

    files = ["riptide"] if riptide_only else ["gbs", "riptide"]
    dfs = [
        pd.read_csv(f"{base_path}/impute_{name}_svm_output_outcome.csv")
        for name in files
    ]

    sex_qc_df = dfs[0] if len(dfs) == 1 else pd.concat(dfs)

    sex_qc_df = tempdf[['Sample_ID', 'Library_ID', 'rfid', 'Sample_Barcode']].merge(
        sex_qc_df, on=['Sample_ID', 'Library_ID'], how='right'
    )

    sex_qc_df['round_calculated'] = current_round

    sex_qc_df['for_training_threshold'] = sex_qc_df.apply(
        lambda row: riptide_for_training_threshold if row['seq_method'] == 'riptide'
        else gbs_for_training_threshold, axis=1
    )
    sex_qc_df['undetermined_threshold'] = sex_qc_df.apply(
        lambda row: riptide_undetermined_threshold if row['seq_method'] == 'riptide'
        else gbs_undetermined_threshold, axis=1
    )

    sex_qc_df = sex_qc_df.rename(columns={
        'Library_ID': 'library_name',
        'Sample_Barcode': 'barcode',
        'sex': 'recorded_sex',
        'chrX': 'chrX_mapped_read_count',
        'chrY': 'chrY_mapped_read_count',
        'allchrom': 'allchrom_mapped_read_count',
        'chrX_perc': 'chrX_mapped_read_percentage_total',
        'chrY_perc': 'chrY_mapped_read_percentage_total',
        'F_proba': 'female_probability',
        'M_proba': 'male_probability',
        'sexqc_outcome': 'sex_qc_outcome'
    })

    sex_qc_df['sex_qc_status'] = sex_qc_df.apply(
        lambda row: row['sex_qc_outcome'] in ['kept_original', 'imputed']
        if pd.notnull(row['sex_qc_outcome']) else None, axis=1
    )

    sex_qc_df['notes'] = None

    sex_qc_df = sex_qc_df[[
        'library_name', 'rfid', 'barcode', 'round_calculated',
        'recorded_sex',
        'chrX_mapped_read_count', 'chrY_mapped_read_count', 'allchrom_mapped_read_count',
        'chrX_mapped_read_percentage_total', 'chrY_mapped_read_percentage_total',
        'female_probability', 'male_probability',
        'for_training_threshold', 'undetermined_threshold',
        'predicted_sex', 'sex_qc_outcome', 'sex_qc_status',
        'notes'
    ]]

    sex_qc_df.to_csv(output_filepath, index=False)
    print(f'Results written out to: {output_filepath}')

    return sex_qc_df


def create_het_qc_log(multi_flowcell_filepath, round_directory, current_round,
                      miss_threshold, gbs_sd_threshold, riptide_sd_threshold,
                      uniqueseqmethod_sd_threshold, output_filepath):
    """
    Purpose: Creates a database-ready table that summarises the heterozygosity QC results.
    """
    tempdf = pd.read_csv(multi_flowcell_filepath, dtype={'rfid': 'str'})
    tempdf = tempdf[['Sample_ID', 'Library_ID', 'rfid', 'Sample_Barcode', 'seq_method']]

    het_dict = {}
    file_prefixes = {
        "raw": f"06_snp_filter/{{chromvar}}_allsexes/{current_round}_{{chromvar}}_allsexes_raw",
        "snp_filtered": f"06_snp_filter/{{chromvar}}_allsexes/{current_round}_{{chromvar}}_allsexes_snp_filtered",
        "snp_and_sample_filtered": f"07_sample_filter/{{chromvar}}_allsexes/{current_round}_{{chromvar}}_allsexes_snp_and_sample_filtered",
    }
    chrom_list = [f"chr{c}" for c in range(1, 21)]

    for key, template in file_prefixes.items():
        for chromvar in chrom_list:
            prefix = template.format(chromvar=chromvar)
            base = f"{round_directory}/output/{prefix}"
            missdf = pd.read_csv(f"{base}.smiss", sep="\t")
            missdf = missdf.rename(columns={'OBS_CT': 'OBS_CT_MISS'})
            hetdf = pd.read_csv(f"{base}.het", sep="\t")
            combdf = hetdf.merge(missdf[["#IID", "MISSING_CT", "OBS_CT_MISS"]], on="#IID", how="inner")
            het_dict[(key, chromvar)] = combdf

    het_qc_df = pd.concat(het_dict).reset_index(names=['filtering_stage', 'chrom', 'old_index']).drop(columns=['old_index'])
    het_qc_df = het_qc_df.merge(tempdf, left_on='#IID', right_on='Sample_ID', how='left')
    het_qc_df = het_qc_df.groupby(['Library_ID', 'Sample_Barcode', 'rfid', '#IID', 'filtering_stage', 'seq_method'])[
        ['O(HOM)', 'E(HOM)', 'OBS_CT', "MISSING_CT", "OBS_CT_MISS"]].sum().reset_index()

    het_qc_df['heterozygosity'] = (het_qc_df['OBS_CT'] - het_qc_df['O(HOM)']) / het_qc_df['OBS_CT']
    het_qc_df['miss_perc'] = (het_qc_df['MISSING_CT'] / het_qc_df['OBS_CT_MISS']) * 100
    het_qc_df['threshold'] = het_qc_df['seq_method'].map({
        'gbs': gbs_sd_threshold, 'riptide': riptide_sd_threshold, 'LSW_lcwgs': uniqueseqmethod_sd_threshold
    })

    mean_dict = het_qc_df.groupby(['filtering_stage', 'seq_method']).heterozygosity.mean().to_dict()
    std_dict = het_qc_df.groupby(['filtering_stage', 'seq_method']).heterozygosity.std().to_dict()

    het_qc_df['key'] = list(zip(het_qc_df['filtering_stage'], het_qc_df['seq_method']))
    het_qc_df['mean'] = het_qc_df['key'].map(mean_dict)
    het_qc_df['std'] = het_qc_df['key'].map(std_dict)
    het_qc_df['lower_limit'] = het_qc_df.apply(lambda row: row['mean'] - (row['std'] * row['threshold']), axis=1)
    het_qc_df['upper_limit'] = het_qc_df.apply(lambda row: row['mean'] + (row['std'] * row['threshold']), axis=1)
    het_qc_df['het_qc_outlier'] = het_qc_df.apply(
        lambda row: row['lower_limit'] <= row['heterozygosity'] <= row['upper_limit'], axis=1
    )
    het_qc_df['miss_qc_outlier'] = het_qc_df.apply(lambda row: row['miss_perc'] < miss_threshold, axis=1)

    het_qc_df['het_qc_outcome'] = het_qc_df.apply(
        lambda row: (
            'fail' if not row['miss_qc_outlier'] and not row['het_qc_outlier']
            else 'suspect' if not row['het_qc_outlier']
            else 'pass'
        ), axis=1
    )
    het_qc_df['het_qc_status'] = het_qc_df.apply(
        lambda row: row['het_qc_outcome'] in ['suspect', 'pass']
        if pd.notnull(row['het_qc_outcome']) else None, axis=1
    )

    het_qc_df = het_qc_df.rename(columns={
        'Library_ID': 'library_name',
        'Sample_Barcode': 'barcode',
        'O(HOM)': 'ohom',
        'E(HOM)': 'ehom',
        'OBS_CT': 'obsct',
        'heterozygosity': 'heterozygosity_proportion',
        'threshold': 'standard_deviations_threshold',
        'het_qc_outcome': 'heterozygosity_qc_outcome',
        'het_qc_status': 'heterozygosity_qc_status',
    })

    het_qc_df['notes'] = None
    het_qc_df['round_calculated'] = current_round

    het_qc_df = het_qc_df[[
        'library_name', 'rfid', 'barcode',
        'round_calculated', 'filtering_stage',
        'ohom', 'ehom', 'obsct',
        'heterozygosity_proportion', 'standard_deviations_threshold',
        'lower_limit', 'upper_limit',
        'heterozygosity_qc_outcome', 'heterozygosity_qc_status',
        'notes'
    ]]

    het_qc_df.to_csv(output_filepath, index=False)
    print(f'Results written out to: {output_filepath}')

    return het_qc_df


def create_missingness_qc_log(multi_flowcell_filepath, round_directory, current_round, miss_threshold, output_filepath):
    """
    Purpose: Creates a database-ready table that summarises the sample missingness QC results.
    """
    tempdf = pd.read_csv(multi_flowcell_filepath, dtype={'rfid': 'str'})
    tempdf = tempdf[['Sample_ID', 'Library_ID', 'rfid', 'Sample_Barcode', 'seq_method']]

    miss_dict = {}
    file_prefixes = {
        "raw": f"06_snp_filter/{{chromvar}}_allsexes/{current_round}_{{chromvar}}_allsexes_raw",
        "snp_filtered": f"06_snp_filter/{{chromvar}}_allsexes/{current_round}_{{chromvar}}_allsexes_snp_filtered",
        "snp_and_sample_filtered": f"07_sample_filter/{{chromvar}}_allsexes/{current_round}_{{chromvar}}_allsexes_snp_and_sample_filtered",
    }
    chrom_list = [f"chr{c}" for c in range(1, 21)] + ["chrX", "chrY", "chrM"]

    for key, template in file_prefixes.items():
        for chromvar in chrom_list:
            prefix = template.format(chromvar=chromvar)
            base = f"{round_directory}/output/{prefix}"
            missdf = pd.read_csv(f"{base}.smiss", sep="\t")
            miss_dict[(key, chromvar)] = missdf

    miss_qc_df = pd.concat(miss_dict).reset_index(names=['filtering_stage', 'chrom', 'old_index']).drop(columns=['old_index'])
    miss_qc_df = miss_qc_df.merge(tempdf, left_on='#IID', right_on='Sample_ID', how='left')
    miss_qc_df = miss_qc_df.groupby(['Library_ID', 'Sample_Barcode', 'rfid', '#IID', 'filtering_stage'])[
        ["MISSING_CT", "OBS_CT"]].sum().reset_index()

    miss_qc_df['miss_perc'] = (miss_qc_df['MISSING_CT'] / miss_qc_df['OBS_CT']) * 100
    miss_qc_df['missingness_threshold'] = miss_threshold
    miss_qc_df['missingness_qc_status'] = miss_qc_df.apply(lambda row: row['miss_perc'] <= miss_threshold, axis=1)

    miss_qc_df.columns = miss_qc_df.columns.str.lower()
    miss_qc_df = miss_qc_df.rename(columns={
        'library_id': 'library_name',
        'sample_barcode': 'barcode',
        'miss_perc': 'missingness_percentage'
    })

    miss_qc_df['notes'] = None
    miss_qc_df['round_calculated'] = current_round

    miss_qc_df = miss_qc_df[[
        'library_name', 'rfid', 'barcode',
        'round_calculated', 'filtering_stage',
        'missing_ct', 'obs_ct',
        'missingness_percentage',
        'missingness_threshold',
        'missingness_qc_status',
        'notes'
    ]]

    miss_qc_df.to_csv(output_filepath, index=False)
    print(f'Results written out to: {output_filepath}')

    return miss_qc_df

##############################################################################################################################################################################

def main():

    print(f"Running Step 8: Generating database QC logs for {current_round}....")

    os.makedirs(output_directory, exist_ok=True)

    run_step(
        'Create drop library log',
        lambda: create_drop_library_log(
            code=code,
            current_round=current_round,
            output_filepath=f"{output_prefix}_drop_library_log.csv"
        )
    )

    run_step(
        'Create drop project log',
        lambda: create_drop_project_log(
            code=code,
            current_round=current_round,
            output_filepath=f"{output_prefix}_drop_project_log.csv"
        )
    )

    run_step(
        'Create metadata QC log',
        lambda: create_metadata_qc(
            sequencing_data=sequencing_data,
            multi_flowcell_filepath=clean_metadata_filepath,
            current_round=current_round,
            output_filepath=f"{output_prefix}_metadata_qc_log.csv"
        )
    )

    run_step(
        'Create demultiplexing QC log',
        lambda: create_demux_qc(
            multi_flowcell_filepath=clean_metadata_filepath,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            current_round=current_round,
            output_filepath=f"{output_prefix}_demux_qc_log.csv"
        )
    )

    run_step(
        'Create BAM QC log',
        lambda: create_bam_qc(
            multi_flowcell_filepath=clean_metadata_filepath,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            genome_size=genome_size,
            current_round=current_round,
            output_filepath=f"{output_prefix}_bam_qc_log.csv"
        )
    )

    run_step(
        'Create sex QC log',
        lambda: create_sex_qc_log(
            multi_flowcell_filepath=clean_metadata_filepath,
            round_directory=round_directory,
            riptide_only=True,
            current_round=current_round,
            riptide_for_training_threshold=0.8,
            riptide_undetermined_threshold=0.9,
            gbs_for_training_threshold=0.98,
            gbs_undetermined_threshold=0.98,
            output_filepath=f"{output_prefix}_original_sex_qc_log.csv"
        )
    )

    run_step(
        'Create heterozygosity QC log',
        lambda: create_het_qc_log(
            multi_flowcell_filepath=clean_metadata_filepath,
            round_directory=round_directory,
            current_round=current_round,
            miss_threshold=10,
            gbs_sd_threshold=4,
            riptide_sd_threshold=5,
            uniqueseqmethod_sd_threshold=5,
            output_filepath=f"{output_prefix}_heterozygosity_qc_log.csv"
        )
    )

    run_step(
        'Create missingness QC log',
        lambda: create_missingness_qc_log(
            multi_flowcell_filepath=clean_metadata_filepath,
            round_directory=round_directory,
            current_round=current_round,
            miss_threshold=10,
            output_filepath=f"{output_prefix}_missingness_qc_log.csv"
        )
    )

    print("\nPipeline complete.")

##############################################################################################################################################################################

if __name__ == "__main__":
    main()
