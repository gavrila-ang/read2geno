import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from utils.validate_paths import validate_paths
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



def create_demographic_table(metadata_filepath, file_prefix):
    
    # Read in metadata
    mdf = pd.read_csv(metadata_filepath)
    
    # List of demographic stats of interest
    templist = ['sex', 'seq_method', 'tissue_type']

    # Dictionary to store wrangled demo stat dfs
    demodict = {}
    
    # Calculate categorical proportions for each demographic metric by project
    for demostat in templist:
    
        mdf[demostat] = mdf[demostat].fillna('missing')
        
        tempdf = mdf.groupby(['Sample_Project'])[demostat]\
            .value_counts(normalize=True, dropna=False)\
            .unstack('Sample_Project')\
            .fillna(0) * 100
    
        # all project calculations
        all_project_series = mdf[demostat].value_counts(normalize=True).fillna(0) * 100
        tempdf2 = pd.DataFrame({'all_projects': all_project_series})
        
        # Store calculations in dict
        demodict[demostat] = pd.concat([tempdf2, tempdf], axis=1, ignore_index=False)
    
    # Create counts of rfid and Sample IDs by project
    tempdf = mdf.groupby(['Sample_Project'])[['rfid', 'Sample_ID']].nunique().T
    # Create counts of rfid and Sample IDs across all projects
    all_project_series = mdf[['rfid', 'Sample_ID']].nunique().T
    tempdf2 = pd.DataFrame({'all_projects': all_project_series})
        # Store calculations in dict
    demodict['counts'] = pd.concat([tempdf2, tempdf], axis=1, ignore_index=False)
    
    # Combine the dataframes stored in the dictionary
        # Create whole numbers for easy viewing
    demographic_table = round(pd.concat(demodict), 1)
    
    # Export the dataframe; keep index columns containing demographic stats info
    demographic_table.to_csv(f"{file_prefix}_demographic_plot.csv", index=True) 


def create_demographic_image(current_round, file_prefix):
        
    tempdf = pd.read_csv(f"{file_prefix}_demographic_plot.csv", index_col=['Unnamed: 0', 'Unnamed: 1']).reset_index(names=['cat', 'cat_var'])
    tempdf = pd.melt(tempdf, id_vars=['cat', 'cat_var'])
    
    groups = ['counts', 'sex', 'seq_method', 'tissue_type']
    
    fig = make_subplots(
        # subplot_titles=groups,
        rows=len(groups),
        cols=1,
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
                text=d["value"].values,#.reshape(len(d["cat_var"].unique()), -1),
                texttemplate="%{text}",
                zmin=0, zmax=(5_000 if i == 1 else 100),                
                colorscale="Greys",
                showscale=(i in [1, 4]),   # single colorbar
                colorbar=dict(
                    title=('count' if i == 1 else 'percentage'),
                    x=1.0,
                    y=(0.425 if i == 4 else 0.95),
                    len=(0.2 if i == 1 else 0.9)
                )
            ),
            row=i,
            col=1
        )
    
    
    fig.update_layout(
        title=f'Figure 1: {current_round} Sample Demographic',
        height=750,
        width=1200,
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
    
    # Export the dataframe
    fig.write_html(f"{file_prefix}_demographic_plot.html")

    print(f"Finished writing image to {file_prefix}_demographic_plot.html")