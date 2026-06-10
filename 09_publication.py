#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 9: Generate publication-ready QC figures
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Cleaned metadata file
    - Database QC log files (from step 8)
    - PLINK QC files (.info, .vmiss, .afreq, .hardy, .eigenvec)
    - Reference panel legend files

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Figure 1: Sample demographic heatmap (.html)
    - Figure 2: Alignment statistics boxplot (.html)
    - Figure 3: Sex QC scatterplots (.html)
    - Figure 4: Sample heterozygosity and missingness scatterplot (.html)
    - Figure 5: Sample QC drop reasons barplot (.html)
    - Figure 6: Per-chromosome PCA scatterplot (.html)
    - Figure 7: SNP statistics heatmaps (.html)
    - Figure 8: SNP filtering outcomes barplot (.html)

Assumptions:
    - Step 8 (08_schema.py) has been run and database logs exist
    - Custom input file directory structure

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
import gc
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from tqdm import tqdm
from time import time

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)
exclude_sampleids_list = user.user_args.read_exclude_sampleids()

# === Load system arguments ===
database_logs_directory = str(sys.argv[1])
report_output_directory = str(sys.argv[2])

print(f"database_logs_directory: {database_logs_directory}")
print(f"report_output_directory: {report_output_directory}")

# === Declare input/output filepaths ===
clean_metadata_filepath = f"{round_directory}/output/01_metqc/flowcell_sheets/{current_round}_metadata_cleaned.csv"
file_prefix = f"{report_output_directory}/{current_round}"
database_prefix = f"{database_logs_directory}/{current_round}"

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

# === Figure 1: Demographic ===

def create_demographic_table(metadata_filepath, file_prefix):
    """
    Purpose: Computes categorical proportions and sample counts per project for the demographic table.
    """
    mdf = pd.read_csv(metadata_filepath)
    templist = ['sex', 'seq_method', 'tissue_type']
    demodict = {}

    for demostat in templist:
        mdf[demostat] = mdf[demostat].fillna('missing')
        tempdf = mdf.groupby(['Sample_Project'])[demostat]\
            .value_counts(normalize=True, dropna=False)\
            .unstack('Sample_Project')\
            .fillna(0) * 100
        all_project_series = mdf[demostat].value_counts(normalize=True).fillna(0) * 100
        tempdf2 = pd.DataFrame({'all_projects': all_project_series})
        demodict[demostat] = pd.concat([tempdf2, tempdf], axis=1, ignore_index=False)

    tempdf = mdf.groupby(['Sample_Project'])[['rfid', 'Sample_ID']].nunique().T
    all_project_series = mdf[['rfid', 'Sample_ID']].nunique().T
    tempdf2 = pd.DataFrame({'all_projects': all_project_series})
    demodict['counts'] = pd.concat([tempdf2, tempdf], axis=1, ignore_index=False)

    demographic_table = round(pd.concat(demodict), 1)
    demographic_table.to_csv(f"{file_prefix}_demographic_plot.csv", index=True)


def create_demographic_image(current_round, file_prefix):
    """
    Purpose: Plots a heatmap of demographic proportions by project (Figure 1).
    """
    tempdf = pd.read_csv(f"{file_prefix}_demographic_plot.csv", index_col=['Unnamed: 0', 'Unnamed: 1']).reset_index(names=['cat', 'cat_var'])
    tempdf = pd.melt(tempdf, id_vars=['cat', 'cat_var'])

    groups = ['counts', 'sex', 'seq_method', 'tissue_type']

    fig = make_subplots(
        rows=len(groups), cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        row_heights=[2, 3, 4, 6]
    )

    for i, g in enumerate(groups, start=1):
        d = tempdf[tempdf["cat"] == g]
        fig.layout[f'yaxis{i}'].title.text = g

        fig.add_trace(
            go.Heatmap(
                x=d["variable"],
                y=d["cat_var"],
                z=d["value"],
                text=d["value"].values,
                texttemplate="%{text}",
                zmin=0, zmax=(5_000 if i == 1 else 100),
                colorscale="Greys",
                showscale=(i in [1, 4]),
                colorbar=dict(
                    title=('count' if i == 1 else 'percentage'),
                    x=1.0,
                    y=(0.425 if i == 4 else 0.95),
                    len=(0.2 if i == 1 else 0.9)
                )
            ),
            row=i, col=1
        )

    fig.update_layout(
        title=f'Figure 1: {current_round} Sample Demographic',
        height=750, width=1200,
        margin=dict(l=120, r=120),
        template='simple_white',
    )

    fig.update_xaxes(ticks="", row=1, col=1)
    fig.update_xaxes(ticks="", row=2, col=1)
    fig.update_xaxes(ticks="", row=3, col=1)
    fig.update_yaxes(title_text="(A) Sample<br>Identifier", row=1, col=1)
    fig.update_yaxes(title_text="(B) Reported<br>Sex", row=2, col=1)
    fig.update_yaxes(title_text="(C) Sequencing<br>Method", row=3, col=1)
    fig.update_yaxes(title_text="(D) Tissue<br>Origin", row=4, col=1)

    fig.write_html(f"{file_prefix}_demographic_plot.html")
    print(f"Finished writing image to {file_prefix}_demographic_plot.html")


##############################################################################################################################################################################

# === Figure 2: Alignment ===

def create_alignment_table(metadata_filepath, bam_qc_log, file_prefix):
    """
    Purpose: Wrangles the BAM QC log into a long-format table suitable for plotting (Figure 2).
    """
    mdf = pd.read_csv(metadata_filepath)
    bdf = pd.read_csv(bam_qc_log)

    mgdf = bdf.merge(
        mdf, left_on=['library_name', 'rfid', 'barcode'],
        right_on=['Library_ID', 'rfid', 'Sample_Barcode'], how='left'
    )

    alignment_table = pd.melt(
        mgdf,
        id_vars=['Sample_ID', 'Library_ID', 'rfid', 'Sample_Project', 'flowcell_id'],
        value_vars=[
            'primary_read_count', 'primary_mapped_perc_primary', 'primary_duplicates_perc_primary',
            'with_itself_and_mate_mapped_perc_primary', 'with_mate_mapped_diff_chr_mapq5_perc_primary',
            'singletons_count', 'duplicate_adjusted_coverage'
        ]
    )

    alignment_table['variable'] = alignment_table['variable'].replace({
        'primary_read_count': 'primary_n',
        'primary_duplicates_perc_primary': 'primary_duplicate_%',
        'primary_mapped_perc_primary': 'primary_mapped_%',
        'with_itself_and_mate_mapped_perc_primary': 'matemapped_%',
        'with_mate_mapped_diff_chr_mapq5_perc_primary': 'matemapped_q5_diffchr_%',
        'singletons_count': 'singletons_n'
    })

    alignment_table['date'] = pd.to_datetime(
        alignment_table['flowcell_id'].str.split("_", expand=True)[0], yearfirst=True, errors='coerce'
    )
    alignment_table['date'] = alignment_table['date'].fillna(pd.to_datetime("2018-01-01", yearfirst=True))

    categorical_order = ['primary_n', 'primary_mapped_%', 'matemapped_%', 'matemapped_q5_diffchr_%',
                         'primary_duplicate_%', 'singletons_n', 'coverage', 'duplicate_adjusted_coverage',
                         'total_mapped_base_count', 'total_mapped_nonduplicate_base_count']
    alignment_table['variable'] = pd.Categorical(alignment_table['variable'], categories=categorical_order, ordered=True)
    alignment_table = alignment_table.sort_values(by=['date', 'variable'])
    alignment_table = alignment_table.query('flowcell_id != "GBS"').reset_index(drop=True)

    alignment_table.to_csv(f"{file_prefix}_alignment_plot.csv", index=False)


def create_alignment_image(file_prefix, current_round, grouping='flowcell_id'):
    """
    Purpose: Plots alignment statistics as faceted boxplots (Figure 2).
    """
    alignment_table = pd.read_csv(f"{file_prefix}_alignment_plot.csv")

    categorical_order = ['primary_n', 'primary_mapped_%', 'primary_duplicate_%', 'duplicate_adjusted_coverage',
                         'matemapped_%', 'matemapped_q5_diffchr_%', 'singletons_n']
    alignment_table['variable'] = pd.Categorical(alignment_table['variable'], categories=categorical_order, ordered=True)
    alignment_table = alignment_table.sort_values(by=['date', 'variable'])

    alignment_table['variable'] = alignment_table['variable'].replace({
        'primary_n': '(A) Primary Read (n)',
        'primary_mapped_%': '(B) Primary Reads Mapped (%)',
        'primary_duplicate_%': '(C) Primary Read Duplicates (%)',
        'duplicate_adjusted_coverage': '(D) Adjusted Coverage (X)',
        'matemapped_%': '(E) Mate Mapped (%)',
        'matemapped_q5_diffchr_%': '(F) Mate Mapped to Diff Chrom Q5 (%)',
        'singletons_n': '(G) Singletons (n)'
    })

    fig = px.box(
        alignment_table, y=grouping, x='value',
        color=grouping,
        facet_col='variable',
        title=f'Figure 2: {current_round} {grouping} Riptide Alignment Statistics',
        height=alignment_table.flowcell_id.nunique() * 25,
        hover_data=['Sample_ID', 'Sample_Project']
    )

    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1].strip()))

    for col_idx, ann in enumerate(fig.layout.annotations, start=1):
        col_name = ann['text']
        col_values = alignment_table.query(f'variable=="{col_name}"').value
        mean_val = col_values.mean()
        std_val = col_values.std()
        upper_3sd = mean_val + (5 * std_val)
        lower_3sd = mean_val - (5 * std_val)
        print(col_name, col_values.median())

        if (col_values < lower_3sd).any().item():
            fig.add_vline(x=lower_3sd, line_width=2, line_color="red", line_dash='dash', opacity=1,
                          annotation_text="lower outlier boundary", annotation_position="top left",
                          annotation_textangle=90, row=1, col=col_idx)

        if (col_values > upper_3sd).any().item():
            fig.add_vline(x=upper_3sd, line_width=2, line_color="red", line_dash='dash', opacity=1,
                          annotation_text="upper outlier boundary", annotation_position="top right",
                          annotation_textangle=90, row=1, col=col_idx)

    fig.update_layout(showlegend=False, template='simple_white', width=300 * len(categorical_order))
    fig.update_xaxes(title=None, matches=None)
    fig.update_yaxes(title='Flowcell ID')
    for col_i in range(2, 8):
        fig.update_yaxes(title="", ticks="", row=1, col=col_i)

    fig.write_html(f"{file_prefix}_{grouping}_alignment_plot.html")
    print(f"Finished writing image to {file_prefix}_{grouping}_alignment_plot.html")


##############################################################################################################################################################################

# === Figure 3: Sex QC ===

def create_sexqc_image(metadata_filepath, sex_qc_log, file_prefix, current_round):
    """
    Purpose: Plots sex QC scatterplots and summary heatmaps (Figure 3).
    """
    sexdf = pd.read_csv(sex_qc_log)
    md = pd.read_csv(metadata_filepath)
    sexdf = sexdf.merge(
        md, left_on=['library_name', 'rfid', 'barcode'],
        right_on=['Library_ID', 'rfid', 'Sample_Barcode'], how='left'
    )

    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.3, 0.7],
        row_heights=[0.5, 0.5],
        horizontal_spacing=0.05,
        vertical_spacing=0.2,
        subplot_titles=(
            "3A: Q20 Alignments to Sex Chromosomes",
            "3C: Sample Validation Summary by Research Projects",
            "3B: Female Class Probability of a Sample",
            "3D: Sample Validation Summary by Select Libraries"
        )
    )

    # --- 3A: Chromosome Q20 alignments ---
    color_dict1 = {'missing': 'black', 'F': 'gold', 'M': 'midnightblue'}
    fig.add_trace(
        go.Scatter(
            x=sexdf['chrX_mapped_read_percentage_total'],
            y=sexdf['chrY_mapped_read_percentage_total'],
            marker_color=sexdf['recorded_sex'].map(color_dict1),
            mode='markers',
            customdata=np.stack([sexdf["Sample_ID"], sexdf["Sample_Project"], sexdf["recorded_sex"]], axis=-1),
            hovertemplate=(
                "chrX %: %{x:.2f}<br>chrY %: %{y:.2f}<br>"
                "Predicted Sex: %{customdata[2]}<br>ID: %{customdata[0]}<br>Project: %{customdata[1]}<extra></extra>"
            )
        ),
        row=1, col=1
    )

    fig.add_annotation(x=2, y=2.25, text="Genetic Male Sample Region", showarrow=False,
                       font=dict(size=14, color="midnightblue"), xref="x", yref="y", row=1, col=1)
    fig.add_annotation(x=10, y=0.5, text="Genetic Female Sample Region", showarrow=False,
                       font=dict(size=14, color="gold"), xref="x", yref="y", row=1, col=1)
    fig.add_shape(type="line", x0=0, x1=11.5, y0=0, y1=1.5,
                  line=dict(color="grey", dash="dash", width=2),
                  label=dict(text="Class Boundary", textposition="end", font=dict(color="black", size=12)),
                  opacity=1, row=1, col=1)
    fig.add_annotation(x=0.15, y=0.98, xref="paper", yref="paper",
                       text="<b>Recorded Sex</b><br><span style='color:gold'>■</span> F<br>"
                            "<span style='color:midnightblue'>■</span> M<br><span style='color:black'>■</span> missing",
                       showarrow=False, align="left", bgcolor="rgba(255,255,255,0.8)",
                       bordercolor="black", borderwidth=1)

    # --- 3B: SVM prediction probabilities ---
    color_dict2 = {
        'record_mismatch': 'red', 'undetermined': 'grey',
        'imputed': 'midnightblue', 'kept_original': 'gold',
    }
    fig.add_trace(
        go.Scatter(
            x=sexdf['female_probability'],
            y=sexdf['chrY_mapped_read_percentage_total'],
            marker_color=sexdf['sex_qc_outcome'].map(color_dict2),
            mode='markers',
            customdata=np.stack([sexdf["Sample_ID"], sexdf["Sample_Project"], sexdf["predicted_sex"]], axis=-1),
            hovertemplate=(
                "chrX %: %{x:.2f}<br>chrY %: %{y:.2f}<br>"
                "Predicted Sex: %{customdata[2]}<br>ID: %{customdata[0]}<br>Project: %{customdata[1]}<extra></extra>"
            )
        ),
        row=2, col=1
    )

    fig.add_vline(x=0.95, line_dash="dash", line_color='grey', line_width=2, opacity=1,
                  annotation_text="female_class_threshold", annotation_position="top right",
                  annotation_textangle=90, row=2, col=1)
    fig.add_vline(x=0.05, line_dash="dash", line_color='grey', line_width=2, opacity=1,
                  annotation_text="male_class_threshold", annotation_position="top right",
                  annotation_textangle=90, row=2, col=1)
    fig.add_annotation(x=0.15, y=0.28, xref="paper", yref="paper",
                       text="<b>SexQC Outcome</b><br><span style='color:gold'>■</span> kept_original<br>"
                            "<span style='color:midnightblue'>■</span> imputed<br>"
                            "<span style='color:red'>■</span> record_mismatch<br>"
                            "<span style='color:grey'>■</span> undetermined",
                       showarrow=False, align="left", bgcolor="rgba(255,255,255,0.8)",
                       bordercolor="black", borderwidth=1)

    # --- 3C: QC outcomes by project ---
    tempdf1 = sexdf.groupby('Sample_Project')['sex_qc_outcome'].value_counts(normalize=True, dropna=False).unstack('Sample_Project').fillna(0)
    tempdf = pd.DataFrame({'all_projects': sexdf['sex_qc_outcome'].value_counts(normalize=True, dropna=False)}).fillna(0)
    bar1df = round(pd.concat([tempdf, tempdf1], axis=1, ignore_index=False) * 100, 2)

    fig.add_trace(
        go.Heatmap(z=bar1df, y=bar1df.index, x=bar1df.columns, text=bar1df.values,
                   texttemplate="%{text}", colorscale="Greys", zmin=0, zmax=100,
                   showscale=True, colorbar=dict(title='percentage', x=1.0, y=0.8, len=0.4)),
        row=1, col=2
    )

    # --- 3D: QC outcomes by library ---
    tempdf1 = sexdf.groupby('Library_ID')['sex_qc_outcome'].value_counts(normalize=False, dropna=False).unstack('Library_ID').fillna(0)
    rows = ['record_mismatch', 'undetermined']
    tempdf1 = tempdf1.loc[:, (tempdf1.loc[rows] != 0).any(axis=0)]
    tempdf = pd.DataFrame({'all_projects': sexdf['sex_qc_outcome'].value_counts(normalize=False, dropna=False).fillna(0)})
    bar2df = pd.concat([tempdf, tempdf1], axis=1, ignore_index=False)

    fig.add_trace(
        go.Heatmap(z=bar2df, y=bar2df.index, x=bar2df.columns, text=bar2df.values,
                   texttemplate='%{text}', colorscale="Greys", zmin=0, zmax=96,
                   showscale=True, colorbar=dict(title='count', x=1.0, y=0.2, len=0.4)),
        row=2, col=2
    )

    fig.update_xaxes(title_text="Female probability of sample", row=1, col=1)
    fig.update_yaxes(title_text="% reads mapped to chrY", row=1, col=1)
    fig.update_xaxes(title_text="% reads mapped to chrX", row=2, col=1)
    fig.update_yaxes(title_text="% reads mapped to chrY", row=2, col=1)
    fig.update_layout(
        xaxis2=dict(tickangle=90, tickfont=dict(size=10)),
        xaxis4=dict(tickangle=90, tickfont=dict(size=10)),
        height=1500,
        width=75 * len(bar2df.columns.tolist()),
        title_text=f"Figure 3: {current_round} of Sample Identity Validation Outcomes",
        showlegend=False,
        template="simple_white"
    )

    fig.write_html(f"{file_prefix}_sexqc_plot.html")
    print(f"Finished writing image to {file_prefix}_sexqc_plot.html")


##############################################################################################################################################################################

# === Figures 4-5: Sample stats ===

def safe_read_csv(path, **kwargs):
    return pd.read_csv(path, **kwargs) if os.path.exists(path) else None


def create_sample_stats_table(metadata_filepath, round_directory, current_round, code, file_prefix, database_prefix):
    """
    Purpose: Merges sample-level QC logs into a single table for downstream plotting (Figures 4-5).
    """
    mdf   = safe_read_csv(metadata_filepath)
    htdf  = safe_read_csv(f"{database_prefix}_heterozygosity_qc_log.csv")
    msdf  = safe_read_csv(f"{database_prefix}_missingness_qc_log.csv")
    osxdf = safe_read_csv(f"{database_prefix}_original_sex_qc_log.csv")
    dlibdf = safe_read_csv(f"{code}/user/drop_library_log.csv")
    bmdf  = safe_read_csv(f"{database_prefix}_bam_qc_log.csv")

    if mdf is not None:
        mdf = mdf[['Sample_ID', 'Library_ID', 'rfid', 'flowcell_id', 'Sample_Barcode',
                   'Sample_Project', 'tissue_type', 'seq_method']]\
            .rename(columns={'Library_ID': 'library_name', 'Sample_Barcode': 'barcode'})

    if htdf is not None:
        htdf = htdf.query('filtering_stage == "snp_filtered"')[[
            'library_name', 'rfid', 'barcode',
            'heterozygosity_proportion', 'lower_limit', 'upper_limit',
            'heterozygosity_qc_outcome', 'heterozygosity_qc_status'
        ]]

    if msdf is not None:
        msdf = msdf.query('filtering_stage == "snp_filtered"')[[
            'library_name', 'rfid', 'barcode', 'missingness_percentage', 'missingness_qc_status'
        ]]

    if bmdf is not None:
        bmdf = bmdf[['library_name', 'rfid', 'barcode',
                     'primary_read_count', 'primary_mapped_perc_primary',
                     'duplicate_adjusted_coverage', 'primary_duplicates_perc_primary']]

    if osxdf is not None:
        osxdf = osxdf[['library_name', 'rfid', 'barcode', 'predicted_sex', 'sex_qc_outcome', 'sex_qc_status']]

        if dlibdf is not None:
            sex_bad_libraries = dlibdf.query(
                'round_dropped == @current_round and notes.str.contains("many_sex_mismatches", na=False)'
            ).library_name.tolist()
            osxdf.loc[
                (osxdf.sex_qc_status == True) & (osxdf.library_name.isin(sex_bad_libraries)),
                'new_sex_qc_outcome'
            ] = 'dropped_library'

        osxdf['new_sex_qc_outcome'] = osxdf.get('new_sex_qc_outcome', osxdf['sex_qc_outcome'])
        osxdf['new_sex_qc_status'] = osxdf['new_sex_qc_outcome'].isin(['kept_original', 'imputed'])

    def safe_merge(left, right, on, how='outer'):
        if left is None: return right
        if right is None: return left
        return left.merge(right, on=on, how=how)

    mgdf = mdf.copy() if mdf is not None else None
    mgdf = safe_merge(mgdf, bmdf,  ['library_name', 'rfid', 'barcode'])
    mgdf = safe_merge(mgdf, osxdf, ['library_name', 'rfid', 'barcode'])
    mgdf = safe_merge(mgdf, htdf,  ['library_name', 'rfid', 'barcode'])
    mgdf = safe_merge(mgdf, msdf,  ['library_name', 'rfid', 'barcode'])

    for _, row in tqdm(mgdf.iterrows(), total=len(mgdf)):
        fail_list = []
        if row['new_sex_qc_outcome'] == 'dropped_library':  fail_list.append("libdropSEX")
        if row['new_sex_qc_outcome'] == 'record_mismatch':  fail_list.append("recmSEX")
        if row['new_sex_qc_outcome'] == 'undetermined':     fail_list.append("udmSEX")
        if row['heterozygosity_qc_outcome'] == 'fail':      fail_list.append("outHET")
        if row['heterozygosity_qc_outcome'] == 'suspect':   fail_list.append("susHET")
        if row["missingness_qc_status"] == False:           fail_list.append("hiMISS")
        fail_reason = '_'.join(fail_list) if fail_list else 'pass'
        mgdf.loc[(mgdf.Sample_ID == row['Sample_ID']), 'fail_reason'] = fail_reason

    mgdf['date'] = mgdf.apply(
        lambda row: pd.to_datetime('2019-01-01') if row['flowcell_id'] == 'GBS'
        else pd.to_datetime(row['flowcell_id'].split('_')[0], yearfirst=True, errors='coerce'), axis=1
    )
    mgdf['year'] = mgdf['date'].dt.year
    mgdf['year_month'] = mgdf['date'].dt.strftime("%Y-%m")
    mgdf = mgdf.sort_values(by='year')

    mgdf.to_csv(f"{file_prefix}_sample_stats_plot.csv", index=False)


def create_sample_heterozygosity_plot(file_prefix, current_round):
    """
    Purpose: Plots sample heterozygosity vs missingness scatterplot (Figure 4).
    """
    df = pd.read_csv(f"{file_prefix}_sample_stats_plot.csv")

    df['missingness_qc_status'] = df['missingness_qc_status'].map({True: 'miss_pass', False: 'miss_fail'})
    df['heterozygosity_qc_outcome'] = df['heterozygosity_qc_outcome'].replace(
        {'pass': 'het_pass', 'fail': 'het_fail', 'suspect': 'het_fail'}
    )
    df['miss_het_qc_outcome'] = df['heterozygosity_qc_outcome'] + ';' + df['missingness_qc_status']

    color_mapping = {
        'het_pass;miss_pass': 'green', 'het_pass;miss_fail': 'brown',
        'het_fail;miss_fail': 'red', 'het_fail;miss_pass': 'orange'
    }

    seq_methods = sorted(df['seq_method'].dropna().unique())

    fig = px.scatter(
        df, x='heterozygosity_proportion', y='missingness_percentage',
        facet_col='seq_method', category_orders={'seq_method': seq_methods},
        color='miss_het_qc_outcome', color_discrete_map=color_mapping,
        hover_name='Sample_ID', hover_data=['library_name', 'Sample_Project']
    )

    fig.for_each_annotation(lambda a: a.update(text=a.text.replace("seq_method=", "").upper()))

    for i in range(1, len(seq_methods) + 1):
        fig.add_hline(y=10, line_dash="dot", line_color='brown', row=1, col=i,
                      showlegend=(i == 1), name='miss qc threshold')

    for i, method in enumerate(seq_methods, start=1):
        sub = df[df['seq_method'] == method]
        for colname in ['lower_limit', 'upper_limit']:
            val = sub[colname].dropna().iloc[0] if not sub[colname].dropna().empty else None
            if val is not None:
                fig.add_vline(x=val, line_dash="dot", line_color='red', row=1, col=i,
                              showlegend=(i == 1 and colname == 'lower_limit'), name='het qc thresholds')

    fig.update_layout(
        title='Figure 4: Missingness and Heterozygosity QC Outcomes',
        yaxis_title="Genotype Missingness (%)", legend_title_text="QC Outcome",
        template='presentation', width=500 * len(seq_methods), height=500
    )

    fig.write_html(f"{file_prefix}_sample_heterozygosity_plot.html")
    print(f"Finished writing image to {file_prefix}_sample_heterozygosity_plot.html")


def create_sample_dropreasons_plot(file_prefix, current_round, normalize=False):
    """
    Purpose: Plots sample QC failure reasons by project (Figure 5).
    """
    df = pd.read_csv(f"{file_prefix}_sample_stats_plot.csv", dtype={'year': str})

    df['fail_reason'] = df['fail_reason'].replace({
        'libdropSEX': 'library dropped',
        'recmSEX': 'sex record mismatch',
        'udmSEX': 'sex undetermined',
        'susHET': 'outlying heterozygosity; low missingness',
        'hiMISS': 'high missingness',
        'outHET_hiMISS': 'outlying heterozygosity; high missingness'
    })

    color_mapping = {
        'library dropped': 'green', 'sex record mismatch': 'pink',
        'sex undetermined': 'darkgrey', 'outlying heterozygosity; low missingness': 'orange',
        'high missingness': 'brown', 'outlying heterozygosity; high missingness': 'red'
    }

    agg = (
        df.groupby('Sample_Project')['fail_reason']
        .value_counts(dropna=False, normalize=normalize)
        .rename('value').reset_index()
        .query('fail_reason != "pass"')
    )
    if normalize:
        agg['value'] = agg['value'] * 100

    fig = px.bar(agg, y='Sample_Project', x='value', text_auto=True,
                 color='fail_reason', color_discrete_map=color_mapping)

    fig.update_layout(
        title='Figure 5 Sample QC Failure Summary by Project',
        yaxis_title='Sample Project',
        xaxis_title='Perc' if normalize else 'Sample Count',
        legend_title_text='QC Failure',
        height=df.Sample_Project.nunique() * 25,
        template='simple_white', width=1500
    )

    fig.write_html(f"{file_prefix}_sample_dropreasons_plot.html")
    print(f"Finished writing image to {file_prefix}_sample_dropreasons_plot.html")


##############################################################################################################################################################################

# === Figure 6: PCA ===

def plot_pca_by_chrom(base_path, current_round, metadata_filepath, file_prefix):
    """
    Purpose: Plots PC1 vs PC2 per chromosome, coloured by Sample_Project (Figure 6).
    """
    chrom_list = [f"chr{i}" for i in range(1, 21)] + ["chrX", "chrY", "chrM"]
    dfs = []

    for chrom in chrom_list:
        path = f"{base_path}/{chrom}_allsexes/{current_round}_{chrom}_allsexes_snp_and_sample_filtered.eigenvec"
        try:
            df = pd.read_csv(path, sep=r"\s+")
        except FileNotFoundError:
            continue
        df = df.rename(columns={"#IID": "IID"})
        df["chrom"] = chrom
        dfs.append(df)

    if not dfs:
        raise ValueError("No eigenvec files found.")

    combdf = pd.concat(dfs, ignore_index=True)

    mdf = pd.read_csv(metadata_filepath)
    mdf = mdf.rename(columns={"Sample_ID": "IID", "Sample_Project": "project_name"})[["IID", "project_name"]]
    combdf = combdf.merge(mdf, on="IID", how="left")

    fig = px.scatter(combdf, x="PC1", y="PC2", facet_col="chrom", facet_col_wrap=5,
                     color="project_name", hover_name="IID", opacity=0.7)

    fig.for_each_annotation(lambda a: a.update(text=a.text.replace("chrom=", "")))
    fig.update_layout(title="Figure 6: PCA (PC1 vs PC2) by Chromosome", template="simple_white", height=1200)

    fig.write_html(f"{file_prefix}_sample_pca_plot.html")
    print(f"Finished writing image to {file_prefix}_sample_pca_plot.html")


##############################################################################################################################################################################

# === Figure 7: SNP statistics heatmaps ===

def create_heatmaps(param, usecols, panel, round_directory, current_round, positions_directory,
                    reference_panels_directory, file_prefix, window_size=100_000):
    """
    Purpose: Builds windowed heatmap data files for a given SNP statistic and panel type (Figure 7).
    """
    concat_list = []
    CHROM_ORDER = [f"chr{i}" for i in range(1, 21)] + ["chrX", "chrY", "chrM"]
    CHROM_CAT = pd.CategoricalDtype(categories=CHROM_ORDER, ordered=True)

    for chrom in CHROM_ORDER:
        print(chrom)
        suffix = '.x' if (param == 'hardy' and chrom == 'chrX') else ''
        formpath = (
            f"{round_directory}/output/06_snp_filter/{chrom}_allsexes/"
            f"{current_round}_{chrom}_allsexes_raw.{param}{suffix}"
        )
        if not os.path.exists(formpath):
            continue

        df = pd.read_csv(formpath, sep="\t", usecols=usecols)
        df = df.rename(columns={"#ID": "ID", "[3]CHROM:[4]POS": "ID"})
        df[['CHROM', 'POS']] = df["ID"].str.split(":", expand=True)
        df['POS'] = df['POS'].astype(int)

        panels = pd.read_csv(
            f"{reference_panels_directory}/vcf1_filtered_{chrom}.legend.gz",
            sep='\\s', engine='python', usecols=[1]
        )['position'].tolist()

        df = df.query('POS.isin(@panels)') if panel == "panel" else df.query('~POS.isin(@panels)')

        if param == 'info':
            df[param] = df['[5]INFO_SCORE'].astype(str).str.split(',').apply(
                lambda row: np.mean([float(x) for x in row])
            )
        elif param == 'vmiss':
            df[param] = df['F_MISS'] * 100
        elif param == 'afreq':
            df[param] = np.minimum(df['ALT_FREQS'], 1.0 - df['ALT_FREQS'])
        elif param == 'hardy':
            df[param] = -np.log10(df['P'].replace(0, np.nan))

        df['window'] = (df['POS'] // window_size) * window_size
        statsdf = df.groupby('window')[param].mean().reset_index().sort_values(by='window')
        statsdf['CHROM'] = chrom
        concat_list.append(statsdf)

        del df, statsdf
        gc.collect()

    concatdf = pd.concat(concat_list, ignore_index=True)
    concatdf['CHROM'] = concatdf['CHROM'].astype(CHROM_CAT)
    concatdf = concatdf.sort_values('CHROM')
    heatdf = concatdf.pivot(index='CHROM', columns='window', values=param)
    heatdf.to_csv(f"{file_prefix}_{param}_{panel}.heatmap")
    print(f"Saved to {file_prefix}_{param}_{panel}.heatmap")


def plot_heatmaps(file_prefix, current_round):
    """
    Purpose: Assembles per-parameter, per-panel heatmaps into a single figure (Figure 7).
    """
    fig = make_subplots(
        rows=4, cols=2,
        vertical_spacing=0.02, horizontal_spacing=0.02,
        subplot_titles=['Reference Panel SNPs', 'Non-Reference Panel SNPs'],
        shared_yaxes=True, shared_xaxes=True
    )

    templist = [(x, y) for x in enumerate(['info', 'vmiss', 'afreq', 'hardy'], start=1) for y in ['panel', 'nonpanel']]

    for (row_num, param), panel in templist:
        print(row_num, param, panel)
        heatdf = pd.read_csv(f"{file_prefix}_{param}_{panel}.heatmap", index_col='CHROM')

        fig.add_trace(
            go.Heatmap(
                z=heatdf.values, x=heatdf.columns, y=heatdf.index,
                colorscale='YlOrRd',
                showscale=(panel == 'nonpanel'),
                colorbar=dict(title='', x=1, y=1.125 - (0.25 * row_num), len=0.25)
            ),
            col=(1 if panel == 'panel' else 2), row=row_num
        )

    fig.update_yaxes(dtick=1)
    fig.update_layout(
        yaxis1_title='(A) STITCH INFO Score',
        yaxis3_title='(B) Sample Missingness (%)',
        yaxis5_title='(C) Minor Allele Frequency',
        yaxis7_title='(D) Hardy-Weinberg Equillibrium (-log10 pval)',
        height=1600, template='simple_white',
        title=f"Figure 7: SNP Statistics",
    )

    fig.write_html(f"{file_prefix}_heatmaps.html")
    print(f"Finished writing image to {file_prefix}_heatmaps.html")


##############################################################################################################################################################################

# === Figure 8: SNP drop reasons ===

def create_snp_stats_table(info_threshold, miss_threshold, maf_threshold, hwe_threshold,
                           positions_directory, reference_panels_directory,
                           round_directory, current_round, file_prefix):
    """
    Purpose: Merges per-SNP QC statistics with panel labels into a single table
    for plotting drop reasons (Figure 8). Samples up to 10k SNPs per chromosome
    to keep the output manageable.
    """
    chrom_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','X','Y','M']

    # === Build panel label lookup (sampled) ===
    positions_dict = {}
    for chrom in chrom_list:
        positions = pd.read_csv(f"{positions_directory}/chr{chrom}_STITCH_expanded_pos", sep='\t', usecols=[1], names=['pos'])
        panels = pd.read_csv(f"{reference_panels_directory}/vcf1_filtered_chr{chrom}.legend.gz",
                             sep='\\s', engine='python', usecols=[1]).position.tolist()
        positions['label'] = positions['pos'].isin(panels)
        positions = positions.sample(n=min(10_000, positions.shape[0]), random_state=98)
        positions_dict[f"chr{chrom}"] = positions

    labels = pd.concat(positions_dict).reset_index(names=['chrom', 'old_idx'], level=0)
    labels['ID'] = labels['chrom'] + ':' + labels['pos'].astype(str)
    labels = labels.set_index('ID').sort_index()
    del positions_dict, panels

    # === Load per-chrom SNP stats from 06_snp_filter ===
    info_dict, miss_dict, maf_dict, hwe_dict = {}, {}, {}, {}

    for chrom in chrom_list:
        base = f"{round_directory}/output/06_snp_filter/chr{chrom}_allsexes/{current_round}_chr{chrom}_allsexes_raw"
        chrom_label = f"chr{chrom}"

        # Info scores
        info_path = f"{base}.info"
        if os.path.exists(info_path):
            df = pd.read_csv(info_path, sep='\t', usecols=['[3]CHROM:[4]POS', '[5]INFO_SCORE'])
            df = df.rename(columns={'[3]CHROM:[4]POS': 'ID', '[5]INFO_SCORE': 'info_score'})
            df['info_score'] = df['info_score'].astype(str).str.split(',').apply(
                lambda x: np.mean([float(i) for i in x])
            )
            info_dict[chrom_label] = df.set_index('ID')['info_score']

        # Missingness
        miss_path = f"{base}.vmiss"
        if os.path.exists(miss_path):
            df = pd.read_csv(miss_path, sep='\t', usecols=['ID', 'F_MISS'])
            df['miss'] = df['F_MISS'] * 100
            miss_dict[chrom_label] = df.set_index('ID')['miss']

        # MAF
        maf_path = f"{base}.afreq"
        if os.path.exists(maf_path):
            df = pd.read_csv(maf_path, sep='\t', usecols=['ID', 'ALT_FREQS'])
            df['maf'] = df['ALT_FREQS'].apply(lambda x: x if x <= 0.5 else 1 - x)
            maf_dict[chrom_label] = df.set_index('ID')['maf']

        # HWE — chrX uses .hardy.x
        hwe_suffix = '.hardy.x' if chrom == 'X' else '.hardy'
        hwe_path = f"{base}{hwe_suffix}"
        if os.path.exists(hwe_path):
            df = pd.read_csv(hwe_path, sep='\t', usecols=['ID', 'P'])
            df['hwe'] = -np.log10(df['P'].replace(0, np.nan))
            hwe_dict[chrom_label] = df.set_index('ID')['hwe']

    # === Concatenate per-chrom series ===
    infoscores = pd.concat(info_dict).reset_index(level=0, drop=True).rename('info_score')
    miss       = pd.concat(miss_dict).reset_index(level=0, drop=True).rename('miss')
    maf        = pd.concat(maf_dict).reset_index(level=0, drop=True).rename('maf')
    hwe        = pd.concat(hwe_dict).reset_index(level=0, drop=True).rename('hwe')

    # === Merge all stats onto sampled label positions ===
    plotdf = labels.copy()
    plotdf = plotdf.join(infoscores, how='left')
    plotdf = plotdf.join(miss,       how='left')
    plotdf = plotdf.join(maf,        how='left')
    plotdf = plotdf.join(hwe,        how='left')

    plotdf = plotdf.reset_index(drop=True)
    plotdf = plotdf.sort_values(by=['chrom', 'pos']).reset_index(drop=True)

    # === Apply QC thresholds ===
    # True = passes QC, False = fails QC (appears in fail_summary)
    plotdf['info_qc'] = (plotdf['info_score'] >= info_threshold) | plotdf['info_score'].isna()
    plotdf['miss_qc'] = (plotdf['miss'] < miss_threshold)        | plotdf['miss'].isna()
    plotdf['maf_qc']  = (plotdf['maf'] > maf_threshold)          | plotdf['maf'].isna()
    plotdf['hwe_qc']  = (plotdf['hwe'] < hwe_threshold)          | plotdf['hwe'].isna()

    # === Create fail summary column ===
    qc_cols = ['info_qc', 'miss_qc', 'maf_qc', 'hwe_qc']
    plotdf['fail_summary'] = plotdf[qc_cols].apply(
        lambda row: '_'.join(col for col in qc_cols if not row[col]),
        axis=1
    )

    # === Create drop status column ===
    plotdf['drop_status'] = plotdf['fail_summary'].str.contains('info|miss', regex=True)

    plotdf.to_csv(f"{file_prefix}_snp_stats_plot.csv", index=False)
    print(f"Saved to {file_prefix}_snp_stats_plot.csv")

    

def create_snp_stats_image(file_prefix, current_round):
    """
    Purpose: Plots the proportion of SNPs failing each QC filter, per chromosome (Figure 8).
    """
    plotdf = pd.read_csv(f"{file_prefix}_snp_stats_plot.csv")

    chrom_df = (
        plotdf.groupby(['chrom', 'label'])['fail_summary']
        .value_counts(dropna=False, normalize=True)
        .reset_index(name='proportion')
    )

    all_df = (
        plotdf.groupby(['label'])['fail_summary']
        .value_counts(dropna=False, normalize=True)
        .reset_index(name='proportion')
    )
    all_df['chrom'] = 'All'

    tempdf = pd.concat([chrom_df, all_df], ignore_index=True)
    tempdf['proportion'] = round(tempdf['proportion'] * 100, 1)

    tempdf['chrom'] = tempdf['chrom'].astype(str)
    tempdf.loc[tempdf['chrom'] != 'All', 'chrom'] = 'chr' + tempdf.loc[tempdf['chrom'] != 'All', 'chrom']
    tempdf['chrom'] = tempdf['chrom'].replace({'chr21': 'chrX', 'chr22': 'chrY', 'chr23': 'chrM'})

    tempdf['fail_summary'] = tempdf['fail_summary'].fillna('all_filters_passed')
    tempdf['fail_summary'] = tempdf['fail_summary'].replace({
        'miss_qc': 'hiMISS', 'maf_qc': 'loMAF',
        'info_qc_miss_qc': 'loINFO; hiMISS', 'info_qc_miss_qc_maf_qc': 'loINFO; hiMISS; loMAF',
        'info_qc_maf_qc': 'loINFO; loMAF', 'info_qc': 'loINFO',
        'info_qc_hwe_qc': 'loINFO; hiHWE', 'hwe_qc': 'hiHWE',
        'info_qc_miss_qc_hwe_qc': 'loINFO; hiMISS; loHWE'
    })

    tempdf = tempdf.query('fail_summary != "all_filters_passed"')

    chrom_order = ['All'] + [f"chr{i}" for i in range(1, 21)] + ['chrX', 'chrY', 'chrM']

    fig = px.bar(
        tempdf, x='proportion', y='chrom', color='fail_summary',
        facet_col='label', facet_col_spacing=0.04, text_auto=True,
        category_orders={'chrom': chrom_order, 'label': [True, False]}
    )

    fig.for_each_annotation(
        lambda a: a.update(
            text=a.text.upper()
            .replace("LABEL=FALSE", "(B) Non-Reference Panel SNPs")
            .replace("LABEL=TRUE", "(A) Reference Panel SNPs")
        )
    )

    fig.update_xaxes(matches=None)
    fig.update_layout(xaxis=dict(autorange="reversed"))
    fig.update_layout(
        title='Figure 8: SNP Filtering Outcomes',
        legend_title_text='QC Failure', height=750, width=2000, template='presentation',
    )
    fig.update_xaxes(title=None)
    fig.add_annotation(
        text='% of total SNPs', x=0.5, y=-0.1, xref='paper', yref='paper',
        showarrow=False, xanchor='center', font=dict(size=20)
    )
    fig.update_yaxes(title='', showticklabels=True, side="right", row=1, col=1)

    fig.write_html(f"{file_prefix}_snp_dropreasons_plot.html")
    print(f"Finished writing image to {file_prefix}_snp_dropreasons_plot.html")


def create_snp_stats_image(file_prefix, current_round):

    plotdf = pd.read_csv(f"{file_prefix}_snp_stats_plot.csv")

    chrom_df = (
        plotdf
        .groupby(['chrom', 'label'])['fail_summary']
        .value_counts(dropna=False, normalize=True)
        .reset_index(name='proportion')
    )

    all_df = (
        plotdf
        .groupby(['label'])['fail_summary']
        .value_counts(dropna=False, normalize=True)
        .reset_index(name='proportion')
    )
    all_df['chrom'] = 'All'

    tempdf = pd.concat([chrom_df, all_df], ignore_index=True)
    tempdf['proportion'] = round(tempdf['proportion'] * 100, 1)

    # === Remove this block — chrom is already prefixed with 'chr' ===
    # tempdf['chrom'] = tempdf['chrom'].astype(str)
    # tempdf.loc[tempdf['chrom'] != 'All', 'chrom'] = 'chr' + tempdf.loc[...]
    # tempdf['chrom'] = tempdf['chrom'].replace({...})

    tempdf['fail_summary'] = tempdf['fail_summary'].fillna('all_filters_passed')
    tempdf['fail_summary'] = tempdf['fail_summary'].replace({
        'miss_qc': 'hiMISS', 'maf_qc': 'loMAF',
        'info_qc_miss_qc': 'loINFO; hiMISS', 'info_qc_miss_qc_maf_qc': 'loINFO; hiMISS; loMAF',
        'info_qc_maf_qc': 'loINFO; loMAF', 'info_qc': 'loINFO',
        'info_qc_hwe_qc': 'loINFO; hiHWE', 'hwe_qc': 'hiHWE',
        'info_qc_miss_qc_hwe_qc': 'loINFO; hiMISS; loHWE'
    })

    tempdf = tempdf.query('fail_summary != "all_filters_passed"')

    # === Correct ordering ===
    chrom_order = ['All'] + [f"chr{i}" for i in range(1, 21)] + ['chrX', 'chrY', 'chrM']

    fig = px.bar(
        tempdf, x='proportion', y='chrom', color='fail_summary',
        facet_col='label', facet_col_spacing=0.04, text_auto=True,
        category_orders={'chrom': chrom_order, 'label': [True, False]}
    )

    fig.for_each_annotation(
        lambda a: a.update(
            text=a.text.upper()
            .replace("LABEL=FALSE", "(B) Non-Reference Panel SNPs")
            .replace("LABEL=TRUE", "(A) Reference Panel SNPs")
        )
    )

    fig.update_xaxes(matches=None)
    fig.update_layout(xaxis=dict(autorange="reversed"))
    fig.update_layout(
        title='Figure 8: SNP Filtering Outcomes',
        legend_title_text='QC Failure', height=750, width=2000, template='presentation',
    )
    fig.update_xaxes(title=None)
    fig.add_annotation(
        text='% of total SNPs', x=0.5, y=-0.1, xref='paper', yref='paper',
        showarrow=False, xanchor='center', font=dict(size=20)
    )
    fig.update_yaxes(title='', showticklabels=True, side="right", row=1, col=1)

    fig.write_html(f"{file_prefix}_snp_dropreasons_plot.html")
    print(f"Finished writing image to {file_prefix}_snp_dropreasons_plot.html")

    

##############################################################################################################################################################################

def main():

    print(f"Running Step 9: Generating publication figures for {current_round}....")

    os.makedirs(report_output_directory, exist_ok=True)

    # === Figure 1: Demographic ===
    run_step(
        'Create demographic table',
        lambda: create_demographic_table(
            metadata_filepath=clean_metadata_filepath,
            file_prefix=file_prefix
        )
    )
    run_step(
        'Plot demographic heatmap (Figure 1)',
        lambda: create_demographic_image(
            current_round=current_round,
            file_prefix=file_prefix
        )
    )

    # === Figure 2: Alignment ===
    run_step(
        'Create alignment table',
        lambda: create_alignment_table(
            metadata_filepath=clean_metadata_filepath,
            bam_qc_log=f"{database_prefix}_bam_qc_log.csv",
            file_prefix=file_prefix
        )
    )
    run_step(
        'Plot alignment statistics (Figure 2)',
        lambda: create_alignment_image(
            file_prefix=file_prefix,
            current_round=current_round,
            grouping='flowcell_id'
        )
    )

    # === Figure 3: Sex QC ===
    run_step(
        'Plot sex QC outcomes (Figure 3)',
        lambda: create_sexqc_image(
            metadata_filepath=clean_metadata_filepath,
            sex_qc_log=f"{database_prefix}_original_sex_qc_log.csv",
            file_prefix=file_prefix,
            current_round=current_round
        )
    )

    # === Figures 4-5: Sample stats ===
    run_step(
        'Create sample stats table',
        lambda: create_sample_stats_table(
            metadata_filepath=clean_metadata_filepath,
            round_directory=round_directory,
            current_round=current_round,
            code=code,
            file_prefix=file_prefix,
            database_prefix=database_prefix
        )
    )
    run_step(
        'Plot sample heterozygosity and missingness (Figure 4)',
        lambda: create_sample_heterozygosity_plot(
            file_prefix=file_prefix,
            current_round=current_round
        )
    )
    run_step(
        'Plot sample drop reasons (Figure 5)',
        lambda: create_sample_dropreasons_plot(
            file_prefix=file_prefix,
            current_round=current_round,
            normalize=False
        )
    )

    # === Figure 6: PCA ===
    run_step(
        'Plot per-chromosome PCA (Figure 6)',
        lambda: plot_pca_by_chrom(
            base_path=f"{round_directory}/output/07_sample_filter",
            current_round=current_round,
            metadata_filepath=clean_metadata_filepath,
            file_prefix=file_prefix
        )
    )

    # === Figure 7: SNP statistics heatmaps ===
    for snp_type in ["panel", "nonpanel"]:
        run_step(
            f'Create {snp_type} SNP heatmap data',
            lambda st=snp_type: [
                create_heatmaps(param='info',  usecols=['[3]CHROM:[4]POS', '[5]INFO_SCORE'], panel=st,
                                round_directory=round_directory, current_round=current_round,
                                positions_directory=positions_directory, reference_panels_directory=reference_panels_directory,
                                file_prefix=file_prefix),
                create_heatmaps(param='vmiss', usecols=['ID', 'F_MISS'], panel=st,
                                round_directory=round_directory, current_round=current_round,
                                positions_directory=positions_directory, reference_panels_directory=reference_panels_directory,
                                file_prefix=file_prefix),
                create_heatmaps(param='afreq', usecols=['ID', 'ALT_FREQS'], panel=st,
                                round_directory=round_directory, current_round=current_round,
                                positions_directory=positions_directory, reference_panels_directory=reference_panels_directory,
                                file_prefix=file_prefix),
                create_heatmaps(param='hardy', usecols=['ID', 'P'], panel=st,
                                round_directory=round_directory, current_round=current_round,
                                positions_directory=positions_directory, reference_panels_directory=reference_panels_directory,
                                file_prefix=file_prefix),
            ]
        )

    run_step(
        'Plot SNP statistics heatmaps (Figure 7)',
        lambda: plot_heatmaps(
            file_prefix=file_prefix,
            current_round=current_round
        )
    )

    # === Figure 8: SNP drop reasons ===
    run_step(
        'Create SNP stats table',
        lambda: create_snp_stats_table(
            info_threshold=0.9,
            miss_threshold=10,
            maf_threshold=0,
            hwe_threshold=10,
            positions_directory=positions_directory,
            reference_panels_directory=reference_panels_directory,
            round_directory=round_directory,
            current_round=current_round,
            file_prefix=file_prefix
        )
    )
    
    run_step(
        'Plot SNP drop reasons (Figure 8)',
        lambda: create_snp_stats_image(
            file_prefix=file_prefix,
            current_round=current_round
        )
    )

    print("\nPipeline complete.")

##############################################################################################################################################################################

if __name__ == "__main__":
    main()
