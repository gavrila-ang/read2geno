#!/usr/bin/env python
import pandas as pd
import subprocess
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter("ignore", category=pd.errors.SettingWithCopyWarning)
import plotly.express as px



##############################################################################################################################################################################

def plot_svm_probs(input_directory, dtype_mapping, stage, seq_method, gbs_sexqc_limit, riptide_sexqc_limit, output_directory):
    
    # Read in svm data
    svm_data = pd.read_csv(f"{input_directory}/{stage}_{seq_method}_svm_output_outcome.csv", dtype=dtype_mapping).reset_index(drop=True)
    
    # Color dictionary for plotting
    color_dict = {'kept_original': '#636EFA', 'undetermined': '#FECB52', 'record_mismatch': '#EF553B', 'imputed': '#00CC96'}

    # Create scatter plot for probabilities
    fig = px.scatter(svm_data, 
                     y='chrY_perc',
                     x='F_proba',
                     color='sexqc_outcome',
                     hover_data=["Library_ID"],
                     width=1000,
                     marginal_x="histogram",
                     marginal_y="histogram",
                     height=600,
                     symbol='predicted_sex',
                     color_discrete_map=color_dict,
                    title=f"Model Confidence in Female Classification vs Percentage of Q20 Mapped Read Counts in ChrY <br>Grouped by Concordance between Recorded and Predicted Sex")

    # Shade in where the undetermined samples are
    if (stage == 'clean') and (seq_method == 'riptide'):
         undetermined_cutoff = 0.8
    elif (stage == 'impute') and (seq_method == 'riptide'):
         undetermined_cutoff = riptide_sexqc_limit         
    elif (stage == 'clean') and (seq_method == 'gbs'):
         undetermined_cutoff = 0.98
    elif (stage == 'impute') and (seq_method == 'gbs'):
         undetermined_cutoff = gbs_sexqc_limit
    else:
        print("ERROR: Invalid stage and seq_method combination.")
        sys.exit(1)

    fig.add_vrect(x0=(1-undetermined_cutoff), x1=undetermined_cutoff, line_width=0, fillcolor="grey", opacity=0.3)

    # Export image    
    fig.write_html(f"{output_directory}/plot_svm_probs_{stage}_{seq_method}.html")
    fig.write_image(f"{output_directory}/plot_svm_probs_{stage}_{seq_method}.png")


##############################################################################################################################################################################

def plot_classic_sexqc(input_directory, dtype_mapping, stage, seq_method, output_directory):

    # Read in svm data
    svm_data = pd.read_csv(f"{input_directory}/{stage}_{seq_method}_svm_output_outcome.csv", dtype=dtype_mapping).reset_index(drop=True)

    # Color dictionary for plotting
    color_dict = {'kept_original': '#636EFA', 'undetermined': '#FECB52', 'record_mismatch': '#EF553B', 'imputed': '#00CC96'}
    
    # Create classic scatter plot 
    fig = px.scatter(svm_data, x='chrX_perc', y='chrY_perc',
                     width=800, height=600,
                     color='sexqc_outcome',
                     color_discrete_map=color_dict
                    )

    # Export image    
    fig.write_html(f"{output_directory}/plot_classic_sexqc_{stage}_{seq_method}.html")
    fig.write_image(f"{output_directory}/plot_classic_sexqc_{stage}_{seq_method}.png")


##############################################################################################################################################################################

def plot_svm_summary(input_directory, dtype_mapping, stage, seq_method, output_directory):

    # Read in svm data
    svm_data = pd.read_csv(f"{input_directory}/{stage}_{seq_method}_svm_output_outcome.csv", dtype=dtype_mapping).reset_index(drop=True)

    for feature in ['Library_ID', 'Sample_Project']:
    
        # Determine count and proportion of samples falling in each sexqc_outcome, per Library
        count_df_prop = svm_data.groupby(feature).sexqc_outcome.value_counts(dropna=False, normalize=True).reset_index()
        count_df_count = svm_data.groupby(feature).sexqc_outcome.value_counts(dropna=False, normalize=False).reset_index()
        count_df_prop = count_df_prop.merge(count_df_count, on=[feature, 'sexqc_outcome'], how='outer')
        
        # Make proportions add up to 100, and round out
        count_df_prop['proportion'] = round(count_df_prop['proportion'] * 100, 2)
        
        # # Drop good libraries from visualisation
        # good_categorical_values = count_df_prop.query('sexqc_outcome == "kept_original" and proportion == 100')[feature].tolist()
        # count_df_prop = count_df_prop.query(f'~{feature}.isin(@good_categorical_values)')
        
        # Pivot so each feature has a column per outcome
        pivot_df = count_df_prop.pivot(index=feature, columns='sexqc_outcome', values='proportion').fillna(0)

        # Order by library having the most record_mismatch, then undetermined, then imputed
        setA = set(['record_mismatch', 'undetermined', 'imputed'])
        setB = set(pivot_df.columns)
        if len(setB.difference(setA)) == 0:
            categorical_order = pivot_df.sort_values(by=['record_mismatch', 'undetermined', 'imputed'], ascending=False).index.tolist()
        else:
            categorical_order = count_df_prop[feature].tolist()

        # State bar color of each outcome
        color_dict = {'kept_original': '#636EFA', 'undetermined': '#FECB52', 'record_mismatch': '#EF553B', 'imputed': '#00CC96'}

        # Create figure
        fig = px.bar(count_df_prop, x='proportion', y=feature, 
                     color='sexqc_outcome', color_discrete_map=color_dict, text='count', 
                     category_orders={feature: categorical_order},
                     width=1000, 
                     height=20 * count_df_prop[feature].nunique(), 
                     title=f"Sex QC outcomes of {feature}")
        fig.update_traces(textposition='inside')

        # Export image
        fig.write_html(f"{output_directory}/plot_svm_{feature}_summary_{stage}_{seq_method}.html")
        fig.write_image(f"{output_directory}/plot_svm_{feature}_summary_{stage}_{seq_method}.png")