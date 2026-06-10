#!/usr/bin/env python
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



##############################################################################################################################################################################

def flagstat_scatter(samtools_flagstat_concat_filepath, metadata_filepath, dtype_mapping, output_directory):
    
    # Read in all multiflowcell file
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Read in all flagstat concatenated statistics file
    flagstat_concat_df = pd.read_csv(samtools_flagstat_concat_filepath, index_col=0).reset_index(drop=True)

    # Merge the multiflowcell file and flagstat file
    flagstat_plotting_df = metadata_df.merge(flagstat_concat_df, how='right', on='Sample_ID')
    
    # Create a scatterplot of primary read counts and primary duplicate percentages
    # With histograms in the margin
    fig = px.scatter(flagstat_plotting_df, x="primary_read_count", y="primary_duplicates_perc_primary", color="Library_ID", 
                     marginal_x="histogram", marginal_y="histogram",
                     height=800, width=1200,
                     hover_name="Sample_ID",
                     hover_data=['flowcell_id'],
                     title="Flagstat Metrics: Sample-Level Duplication Rate vs. Primary Read Count, Coloured by Library ID",
                     labels={"primary_read_count": "Primary Reads (Count)", 
                             "primary_duplicates_perc_primary": "Duplicate Reads (% of Primary)"}
              )
    
    # Change axes range and don't show legend (there's too many libraries anyways)
    x_max = flagstat_plotting_df.primary_read_count.max() * 1.1
    y_max = flagstat_plotting_df.primary_duplicates_perc_primary.max() * 1.1
    fig.update_xaxes(range=[0, x_max], row=1, col=1)
    fig.update_yaxes(range=[0, y_max], row=1, col=1)
    fig.update_layout(showlegend=False)

    # Add threshold line for high duplication rates
    hrect_max = flagstat_plotting_df.primary_duplicates_perc_primary.max() * 1.1
    fig.add_hline(y=20, line_color="red", line_width=0.75, line_dash='dash', row=1, col=1)
    fig.add_hline(y=20, line_color="red", line_width=0.75, line_dash='dash', row=1, col=2)
    
    # Add threshold line for minimum read count
    fig.add_vline(x=2_000_000, line_color="red", line_width=0.75, line_dash='dash', row=1, col=1)
    fig.add_vline(x=2_000_000, line_color="red", line_width=0.75, line_dash='dash', row=2, col=1)
    
    # Add rectangle to the plot to highlight samples with high duplication rates and low read counts
    fig.add_shape(
        type="rect",
        x0=0, x1=2000000, y0=20, y1=y_max,
        line=dict(width=0),
        fillcolor="red", opacity=0.3
    )
    
    # Create plotting files
    fig.write_html(f"{output_directory}/samtools_flagstat.html")
    fig.write_image(f"{output_directory}/samtools_flagstat.png")


##############################################################################################################################################################################

def flagstat_boxplot(samtools_flagstat_concat_filepath, metadata_filepath, output_directory, dtype_mapping, order_by="mean"):
    
    flagstat_stats = ["primary_read_count", "primary_perc_total",
                      "secondary_read_count", "secondary_perc_total",
                      "primary_duplicates_read_count", "primary_duplicates_perc_primary",
                      "primary_mapped_read_count", "primary_mapped_perc_primary",
                      "with_itself_and_mate_mapped_read_count", "with_itself_and_mate_mapped_perc_primary",
                      "with_mate_mapped_diff_chr_mapq5_read_count", "with_mate_mapped_diff_chr_mapq5_perc_primary"
                     ]
    
    # Read in all multiflowcell file
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    if order_by == "date":
        # metadata_df['date'] = pd.to_datetime(metadata_df.flowcell_id.str.split("_", expand=True)[0], 
        #                                         yearfirst=True)
        metadata_df['date'] = metadata_df.apply(
            lambda row: pd.to_datetime(row.flowcell_id.split("_")[0], yearfirst=True, errors='coerce')
            if row.flowcell_id != "GBS" else pd.NaT,
            axis=1
        )

    
    # Read in all flagstat concatenated statistics file
    flagstat_concat_df = pd.read_csv(samtools_flagstat_concat_filepath, index_col=0).reset_index(drop=True)

    # Merge the multiflowcell file and flagstat file
    flagstat_plotting_df = metadata_df.merge(flagstat_concat_df, how='right', on='Sample_ID')
    
    # Create figure with number of subplots based on number of flagstat stats of interest
    fig = make_subplots(
        cols=len(flagstat_stats), 
        rows=2,
        row_heights=[0.2 , 0.8],
        vertical_spacing=0.05,
        horizontal_spacing=0.04,
        subplot_titles=[stat for stat in flagstat_stats] + [stat for stat in flagstat_stats],
        shared_xaxes=False
    )

    # Create 2x subplots (boxplot + histogram) for each flagstat
    for i, stat in enumerate(flagstat_stats):
        
        ########################################################################################################################
        
        if order_by == "mean":
            ## Calculate mean per Library_ID for ordering
            order = flagstat_plotting_df.groupby("Library_ID")[stat].mean().sort_values().index.tolist()
        elif order_by == "date":
            order = flagstat_plotting_df[["Library_ID", "date"]].drop_duplicates().sort_values("date")["Library_ID"].tolist()            
        
        ## Set ordered categorical
        ordered_df = flagstat_plotting_df.copy()
        ordered_df["Library_ID"] = pd.Categorical(ordered_df["Library_ID"], categories=order, ordered=True)
        
        ## Create boxplot
        fig_box = px.box(ordered_df, 
                         x=stat, y="Library_ID",
                         color="flowcell_id",
                         hover_name="Sample_ID", hover_data=['flowcell_id'],
                         title="Flagstat Metrics: Sample-Level Duplication Rate vs. Primary Read Count, Coloured by Library ID",
                        )
        
        ## Add each aspect of the boxplot to the main figure
        for trace in fig_box.data:
            fig.add_trace(trace, row=2, col=i+1)
        
        
        ## Order the categorical y-axis of the boxplot
        fig.update_yaxes(categoryorder='array', categoryarray=order,
                         row=2, col=i+1,
                         tickangle=0
                        )

        ## Set the axes range to enable proper rendering of vrect
        x_max = flagstat_plotting_df[stat].max() * 1.1 
        fig.update_xaxes(range=[0, x_max], row=2, col=i+1)
        fig.update_xaxes(range=[0, x_max], row=1, col=i+1)
        fig.update_layout(showlegend=False)
        
        ########################################################################################################################
        
        ## Create histogram        
        fig_hist = px.histogram(flagstat_plotting_df,
                                x=stat,
                                color_discrete_sequence=['gray']
                               )

        ## Add each aspect of the histogram to the main figure
        for trace in fig_hist.data:
            trace.showlegend = False
            fig.add_trace(trace, row=1, col=i+1)

        ########################################################################################################################
            
        ## Draw global (across libraries) lines
        mean_val = flagstat_plotting_df[stat].mean()
        std_val = flagstat_plotting_df[stat].std()
        fig.add_vline(x=mean_val, line_color="black", line_width=1, col=i+1, row=1)
        fig.add_vline(x=mean_val, line_color="black", line_width=1, col=i+1, row=2)

        lower_bound = mean_val - 3 * std_val
        if not lower_bound <= 0:
            fig.add_vline(x=lower_bound, line_color="red", line_dash='dash', line_width=1, col=i+1, row=1)
            fig.add_vline(x=lower_bound, line_color="red", line_dash='dash', line_width=1, col=i+1, row=2)

        upper_bound = mean_val + 3 * std_val  
        fig.add_vline(x=upper_bound, line_color="red", line_dash='dash', line_width=1, col=i+1, row=1)
        fig.add_vline(x=upper_bound, line_color="red", line_dash='dash', line_width=1, col=i+1, row=2)
        
        
        ########################################################################################################################
        
        ## Draw guide rectangles to identify problematic sampels
        if stat in ["secondary_read_count","secondary_perc_total",
                      "primary_duplicates_read_count", "primary_duplicates_perc_primary", 
                      "with_mate_mapped_diff_chr_mapq5_read_count", "with_mate_mapped_diff_chr_mapq5_perc_primary"]:
            fig.add_vrect(x0=upper_bound, x1=x_max, line_width=0, fillcolor="red", opacity=0.3, col=i+1, row=1)
            fig.add_vrect(x0=upper_bound, x1=x_max, line_width=0, fillcolor="red", opacity=0.3, col=i+1, row=2)
        else:
            fig.add_vrect(x0=0, x1=lower_bound, line_width=0, fillcolor="red", opacity=0.3, col=i+1, row=1)
            fig.add_vrect(x0=0, x1=lower_bound, line_width=0, fillcolor="red", opacity=0.3, col=i+1, row=2)


    # Update layout parameters based on number of flagstat stats
    fig.update_layout(
        width=400 * len(flagstat_stats),
        height=50 * flagstat_plotting_df.Library_ID.nunique(),
        showlegend=False,
        title_text="Flagstat Metrics: Mean-Ordered with Side Histograms"
    )

    fig.write_image(f"{output_directory}/flagstat_boxplot_{order_by}_ordered.png")
    fig.write_html(f"{output_directory}/flagstat_boxplot_{order_by}_ordered.html")


##############################################################################################################################################################################

def idxstats_boxplot(samtools_idxstats_concat_filepath, metadata_filepath, output_directory, dtype_mapping):
    
    idxstats_stats = ["count_mapped_read_segments", "count_unmapped_read_segments", "percentage_unmapped_read_segments"]
    
    # Read in all multiflowcell file
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Read in all idxstat concatenated statistics file
    idxstats_concat_df = pd.read_csv(samtools_idxstats_concat_filepath, index_col=0).reset_index(drop=True)
    idxstats_concat_df_melted = pd.melt(idxstats_concat_df, id_vars=["chrom", "idxstat_name"], var_name="Sample_ID", value_name="idxstat_value")
    idxstats_concat_df_pivoted = idxstats_concat_df_melted.pivot(index=["chrom", "Sample_ID"], columns='idxstat_name', values='idxstat_value').reset_index()

    idxstats_concat_df_pivoted["percentage_unmapped_read_segments"] = (idxstats_concat_df_pivoted["count_unmapped_read_segments"] / (idxstats_concat_df_pivoted["count_mapped_read_segments"] + idxstats_concat_df_pivoted["count_unmapped_read_segments"])) * 100
    
    
    # Merge the multiflowcell file and idxstat file
    idxstats_plotting_df = metadata_df.merge(idxstats_concat_df_pivoted, how='right', on='Sample_ID')
    
    # Create figure with number of subplots based on number of flagstat stats of interest
    fig = make_subplots(
        rows=len(idxstats_stats),
        vertical_spacing=0.1,
        subplot_titles=idxstats_stats,
        shared_xaxes=False
    )
    
    # State chromosome order for x-axis 
    order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
            "chrX", "chrY", "chrM"]

    for i, stat in enumerate(idxstats_stats):
        
        ########################################################################################################################

        # Set ordered categorical
        ordered_df = idxstats_plotting_df.copy()
        ordered_df["chrom"] = pd.Categorical(ordered_df["chrom"], categories=order, ordered=True)

        # Create boxplot
        fig_box = px.box(ordered_df,
                         y=stat, x="chrom",
                         color="chrom"
                        )
        for trace in fig_box.data:
            fig.add_trace(trace, row=i+1, col=1)

        # Apply ordering to this subplot's x-axis
        fig.update_xaxes(categoryorder='array', categoryarray=order,
                         row=i+1, col=1
                        )
        
        ## Set the axes range to enable proper rendering of vrect
        y_max = idxstats_plotting_df[stat].max() * 1.1 
        fig.update_yaxes(range=[0, y_max], row=i+1, col=1)
        fig.update_layout(showlegend=False)
        
        ########################################################################################################################

        # Draw global lines
        mean_val = idxstats_plotting_df[stat].mean()
        std_val = idxstats_plotting_df[stat].std()
        fig.add_hline(y=mean_val, line_color="black", line_dash="dash", line_width=1, row=i+1, col=1)

        lower_bound = mean_val - 3 * std_val
        if not lower_bound <= 0:
            fig.add_hline(y=lower_bound, line_color="red", line_width=1, row=i+1, col=1)

        upper_bound = mean_val + 3 * std_val
        fig.add_hline(y=upper_bound, line_color="red", line_width=1, row=i+1, col=1)
        
        ########################################################################################################################

        ## Draw guide rectangles to identify problematic samples
        fig.add_hrect(y0=0, y1=lower_bound, line_width=0, fillcolor="red", opacity=0.3, row=i+1, col=1)

    # Update layout
    fig.update_layout(
        height=300 * len(idxstats_stats),
        width=50*len(order),
        showlegend=False,
        title_text="Idxstats Metrics: Coloured by Library ID",
    )

    fig.write_image(f"{output_directory}/idxstats_histogram.png")
    fig.write_html(f"{output_directory}/idxstats_histogram.html")

##############################################################################################################################################################################

def mosdepth_boxplot(mosdepth_concat_filepath, metadata_filepath, output_directory, dtype_mapping):
    
    mosdepth_stats = ["bases", "length", "max", "mean", "min"]
    
    # Read in all multiflowcell file
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)
    
    # Read in all mosdepth concatenated statistics file
    mosdepth_concat_df = pd.read_csv(mosdepth_concat_filepath, index_col=0).reset_index(drop=True)
    mosdepth_concat_melted = pd.melt(mosdepth_concat_df, id_vars=["chrom", "mosdepth_stat_name"], var_name="Sample_ID", value_name="mosdepth_value")
    mosdepth_concat_pivoted = mosdepth_concat_melted.pivot(index=["chrom", "Sample_ID"], columns='mosdepth_stat_name', values='mosdepth_value').reset_index()
    
    # Merge the multiflowcell file and mosdepth file
    mosdepth_plotting_df = metadata_df.merge(mosdepth_concat_pivoted, how='right', on='Sample_ID')
    
    
    # Create figure with number of subplots based on number of flagstat stats of interest
    fig = make_subplots(
        rows=len(mosdepth_stats),
        vertical_spacing=0.1,
        subplot_titles=mosdepth_stats,
        shared_xaxes=False
    )
    
    # State chromosome order for x-axis 
    order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
            "chrX", "chrY", "chrM"]

    for i, stat in enumerate(mosdepth_stats):
        
        ########################################################################################################################

        # Set ordered categoricals
        ordered_df = mosdepth_plotting_df.copy()
        ordered_df["chrom"] = pd.Categorical(ordered_df["chrom"], categories=order, ordered=True)

        # Create boxplot
        fig_box = px.box(ordered_df,
                         y=stat, x="chrom",
                         color="chrom"
                        )
        for trace in fig_box.data:
            fig.add_trace(trace, row=i+1, col=1)

        # Apply ordering to this subplot's x-axis
        fig.update_xaxes(categoryorder='array', categoryarray=order,
                         row=i+1, col=1
                        )
        
        ## Set the axes range to enable proper rendering of vrect
        y_max = mosdepth_plotting_df[stat].max() * 1.1 
        fig.update_yaxes(range=[0, y_max], row=i+1, col=1)
        fig.update_layout(showlegend=False)
        
        ########################################################################################################################

        # Draw global lines
        mean_val = mosdepth_plotting_df[stat].mean()
        std_val = mosdepth_plotting_df[stat].std()
        fig.add_hline(y=mean_val, line_color="black", line_dash="dash", line_width=1, row=i+1, col=1)

        lower_bound = mean_val - 3 * std_val
        if not lower_bound <= 0:
            fig.add_hline(y=lower_bound, line_color="red", line_width=1, row=i+1, col=1)

        upper_bound = mean_val + 3 * std_val
        fig.add_hline(y=upper_bound, line_color="red", line_width=1, row=i+1, col=1)
        
        ########################################################################################################################

        ## Draw guide rectangles to identify problematic samples
        fig.add_hrect(y0=0, y1=lower_bound, line_width=0, fillcolor="red", opacity=0.3, row=i+1, col=1)

    # Update layout
    fig.update_layout(
        height=300 * len(mosdepth_stats),
        width=50*len(order),
        showlegend=False,
        title_text="Mosdepth Metrics: Coloured by Library ID",
    )

    fig.write_image(f"{output_directory}/mosdepth_histogram.png")
    fig.write_html(f"{output_directory}/mosdepth_histogram.html")


##############################################################################################################################################################################

