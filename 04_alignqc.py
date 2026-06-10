#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 4: Quality control of alignments and sample validation
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Cleaned metadata file
    - Alignment quality control files (.csv, .png, .fastqc, .flagstat, etc.)

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Quality control files (.csv, .png)
    - STITCH input files (sample list, bam list, input file)
    - Cleaned pedigree file

Assumptions:
    - Custom input file directory structure
    - Consistent QC thresholds across all samples

"""

##############################################################################################################################################################################

# === Load python packages ===
import user
import utils
import subprocess
import pandas as pd
import sys
from time import time


# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Declare output directories ===
concat_dir = f"{round_directory}/output/04_align_qc/concat"
html_figures_dir = f"{round_directory}/output/04_align_qc/html_figures"
sexqc_dir = f"{round_directory}/output/04_align_qc/sex_qc"

# === Load system arguments ===
run_type = str(sys.argv[1]) 
gbs_sexqc_limit = float(sys.argv[2]) 
riptide_sexqc_limit = float(sys.argv[3]) 
skip_sexqc = str(sys.argv[4])
stitch_directory = str(sys.argv[5])

print(f"run_type : {run_type}")
print(f"gbs_sexqc_limit : {gbs_sexqc_limit}")
print(f"riptide_sexqc_limit : {riptide_sexqc_limit}")
print(f"skip_sexqc : {skip_sexqc}")

# === Declare input/output filepaths ===
clean_metadata_filepath = f"{round_directory}/output/01_metqc/flowcell_sheets/{current_round}_metadata_cleaned.csv"

# === Declare data type of special columns ===
dtype_mapping={'rfid':str}


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
    print(f"Time elapsed: {elapsed}")

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

def alignqc_block():

    print(f"Running Step 4: Quality control of alignments and sample validation")

    run_step(
        'Aggregates samtools flagstat QC results from multiple samples into a plottable dataframe',
        lambda: utils.align_qc.samtools_flagstat_concat(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            genome_shorthand=genome_shorthand,
            flowcells_directory=flowcells_directory,
            output_filepath=f'{concat_dir}/samtools_flagstat_concat.csv',
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Scatterplot of duplicate read percentages vs primary read counts',
        lambda: utils.align_plot.flagstat_scatter(
            metadata_filepath=clean_metadata_filepath,
            samtools_flagstat_concat_filepath=f"{concat_dir}/samtools_flagstat_concat.csv",
            output_directory=html_figures_dir,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Vertical boxplot of alignment stats grouped by libraries and ordered by mean',
        lambda: utils.align_plot.flagstat_boxplot(
            metadata_filepath=clean_metadata_filepath,
            samtools_flagstat_concat_filepath=f"{concat_dir}/samtools_flagstat_concat.csv",
            order_by='mean',
            output_directory=html_figures_dir,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Vertical boxplot of alignment stats grouped by libraries and ordered by date',
        lambda: utils.align_plot.flagstat_boxplot(
            metadata_filepath=clean_metadata_filepath,
            samtools_flagstat_concat_filepath=f"{concat_dir}/samtools_flagstat_concat.csv",
            order_by='date',
            output_directory=html_figures_dir,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Aggregates samtools idxstats QC results from multiple samples into a plottable dataframe',
        lambda: utils.align_qc.samtools_idxstats_concat(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            genome_shorthand=genome_shorthand,
            flowcells_directory=flowcells_directory,
            output_filepath=f'{concat_dir}/samtools_idxstats_concat.csv',
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Boxplot of mapped base counts by chromosome',
        lambda: utils.align_plot.idxstats_boxplot(
            metadata_filepath=clean_metadata_filepath,
            samtools_idxstats_concat_filepath=f'{concat_dir}/samtools_idxstats_concat.csv',
            output_directory=html_figures_dir,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Aggregates mosdepth QC results from multiple samples into a plottable dataframe',
        lambda: utils.align_qc.mosdepth_concat(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            genome_shorthand=genome_shorthand,
            flowcells_directory=flowcells_directory,
            output_directory=concat_dir,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Boxplot of mapped read counts by chromosome',
        lambda: utils.align_plot.mosdepth_boxplot(
            metadata_filepath=clean_metadata_filepath,
            mosdepth_concat_filepath=f'{concat_dir}/mosdepth_concat_excluded_duplicates.csv',
            output_directory=html_figures_dir,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Aggregates samtools depth QC results from multiple samples into a plottable dataframe',
        lambda: utils.align_qc.samtools_depth_concat(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            genome_shorthand=genome_shorthand,
            flowcells_directory=flowcells_directory,
            output_directory=concat_dir,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Aggregates samtools view results of mapped read counts by chromosome into a plottable dataframe',
        lambda: utils.align_qc.sex_svm_concat(
            metadata_filepath=clean_metadata_filepath,
            exclude_sampleids_list=exclude_sampleids_list,
            genome_shorthand=genome_shorthand,
            flowcells_directory=flowcells_directory,
            output_filepath=f'{concat_dir}/samtools_view_xy_concat.csv',
            dtype_mapping=dtype_mapping            
        )
    )


def sexqc_block(stage, seq_method):

    run_step(
        'Create datasets to run linear SVM model and identify a good training set.',
        lambda: utils.sex_qc.create_svm_input(
            metadata_filepath=clean_metadata_filepath,
            samtools_view_filepath=f"{concat_dir}/samtools_view_xy_concat.csv",
            exclude_sampleids_list=exclude_sampleids_list,
            dtype_mapping=dtype_mapping,
            output_directory=sexqc_dir,
            stage=stage,
            seq_method=seq_method
        )
    )

    run_step(
        'Run linear SVM model.',
        lambda: utils.sex_qc.run_svm_linear(
            input_directory=sexqc_dir,
            dtype_mapping=dtype_mapping,
            output_directory=sexqc_dir,
            downsample='no',
            stage=stage,
            seq_method=seq_method            
        )
    )

    run_step(
        'Determine outcomes of sample validation by comparing predicted and recorded sex values.',
        lambda: utils.sex_qc.create_sexqc_outcome_col(
            input_directory=sexqc_dir,
            dtype_mapping=dtype_mapping,
            gbs_sexqc_limit=gbs_sexqc_limit,
            riptide_sexqc_limit=riptide_sexqc_limit,
            output_directory=sexqc_dir,
            stage=stage,
            seq_method=seq_method            
        )
    )

    run_step(
        'Create a good training set to train a robust linear SVM model.',
        lambda: utils.sex_qc.create_svm_clean_input(
            input_directory=sexqc_dir,
            dtype_mapping=dtype_mapping,
            output_directory=sexqc_dir,
            stage=stage,
            seq_method=seq_method            
        )
    )

    run_step(
        'Scatterplot of probabilities obtained from linear SVM model.',
        lambda: utils.sex_qc_plot.plot_svm_probs(
            input_directory=sexqc_dir,
            dtype_mapping=dtype_mapping,
            gbs_sexqc_limit=gbs_sexqc_limit,
            riptide_sexqc_limit=riptide_sexqc_limit,
            output_directory=sexqc_dir,
            stage=stage,
            seq_method=seq_method            
        )
    )

    run_step(
        'Scatterplot of chrX and chrY mapped read counts.',
        lambda: utils.sex_qc_plot.plot_classic_sexqc(
            input_directory=sexqc_dir,
            dtype_mapping=dtype_mapping,
            output_directory=sexqc_dir,
            stage=stage,
            seq_method=seq_method            
        )
    )

    run_step(
        'Barplot of sample validation outcomes by Library ID and Sample Project.',
        lambda: utils.sex_qc_plot.plot_svm_summary(
            input_directory=sexqc_dir,
            dtype_mapping=dtype_mapping,
            output_directory=sexqc_dir,
            stage=stage,
            seq_method=seq_method            
        )
    )



def stitch_block():

    # run_step(
    #     'Create sample validation summary file, list of samples for STITCH, and list of samples for exclusion.',
    #     lambda: utils.sex_qc.svm_final_outputs(
    #         metadata_filepath=clean_metadata_filepath,
    #         input_directory=sexqc_dir,
    #         dtype_mapping=dtype_mapping,
    #         current_round=current_round,
    #         skip_sexqc=skip_sexqc,
    #         output_directory=sexqc_dir
    #     )
    # )

    run_step(
        'Splits chromosomes into equal-sized chunks for STITCH arrays',
        lambda: utils.stitch_split_chr.run_stitch_split_chr(
            code_path=f"{code}/utils",
            positions_directory=positions_directory,
            output_directory=f"{stitch_directory}/chunk_positions"
        )
    )

    run_step(
        'Creates sex-specific sample lists and bam lists for STITCH runs.',
        lambda: utils.create_input_files.create_stitch_lists(
            metadata_filepath=clean_metadata_filepath,
            flowcells_directory=flowcells_directory,
            genome_shorthand=genome_shorthand,
            current_round=current_round,
            stitch_samples_filepath=f"{round_directory}/output/04_align_qc/sex_qc/{current_round}_stitch_samples.csv",
            output_directory=f"{stitch_directory}/stitch_lists",
            log_directory=f"{round_directory}/logs/04_align_qc",
            exclude_sampleids_list=exclude_sampleids_list,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Creates a STITCH input file containing STITCH parameters based on sex and chromosome.',
        lambda: utils.create_input_files.create_stitch_input(
            chunk_positions_directory=f"{stitch_directory}/chunk_positions",
            stitch_list_directory=f"{stitch_directory}/stitch_lists",
            positions_directory=positions_directory,
            reference_panels_directory=reference_panels_directory,
            current_round=current_round,
            chunk_vcfs_directory=f"{stitch_directory}/stitch_chrom_chunks_vcfs",
            output_directory=f"{stitch_directory}"
        )
    )



def helper_block():

    run_step(
        'Creates a file containing sex data to enable VCF file conversion into PLINK.',
        lambda: utils.create_plink_sex_data(
            stitch_samples_filepath=f"{round_directory}/output/04_align_qc/sex_qc/{current_round}_stitch_samples.csv",
            output_directory=f"{round_directory}/output/09_plink_qc",
            current_round=current_round,
            dtype_mapping=dtype_mapping
        )
    )

    run_step(
        'Creates a file containing pedigree information wrangled from the metadata.',
        lambda: utils.create_pedigree_file(
            metadata_filepath=clean_metadata_filepath,
            imputed_sex_filepath=f"{sexqc_dir}/{current_round}_stitch_samples.csv",
            current_round=current_round,
            exclude_sampleids_list=exclude_sampleids_list,
            output_directory=f"{round_directory}/output/07_snp_filter",
            dtype_mapping=dtype_mapping
        )
    )

    
##############################################################################################################################################################################

def main():

    # === Concatenate alignment quality control files and plot ===
    alignqc_block()

    # === Perform sample validation using observed and genetic sex comparisons ===
    templist = [("clean", "riptide"), ("impute", "riptide"), ("clean", "gbs"), ("impute", "gbs")]
    templist = [("clean", "riptide"), ("impute", "riptide")]
    for stage, seq_method in templist:
        sexqc_block(
            stage=stage, 
            seq_method=seq_method
            )

    # === Create various input files for STITCH ===
    stitch_block()

    # === Create various input files for downstream processes ===
    helper_block()


##############################################################################################################################################################################
        
if __name__ == "__main__":
    main()