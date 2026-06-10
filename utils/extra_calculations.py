# Calculate seq-method specific outlier limits 
def outlier_limits(series, n):
    series_mean = series.mean()
    series_std = series.std()
    
    lower_limit = series_mean - (series_std * n)
    upper_limit = series_mean + (series_std * n)
    
    return lower_limit, upper_limit


# Create a heterogzygosity outlier status function
def het_outlier(het_value, lower_limit, upper_limit, miss_value, miss_limit):

    status_val = None
    
    if (het_value >= upper_limit) and (miss_value >= miss_limit):
        status_val = 'fail'
    elif het_value <= lower_limit and (miss_value >= miss_limit):
        status_val = 'fail'
    elif (het_value >= upper_limit) and (miss_value < miss_limit):
        status_val = 'suspect'
    elif het_value <= lower_limit and (miss_value < miss_limit):
        status_val = 'suspect'        
    else:
        status_val = 'pass'   
        
    return status_val	


def process_heterozygosity_file(prefix):
    
    # Load missingness file
    miss_df = pd.read_csv(f"{prefix}.smiss", sep='\t')
    miss_df = miss_df.rename(columns={'IID': 'Sample_ID'})
    miss_df['missingness_perc'] = miss_df['F_MISS'] * 100
        
    # Load processed heterozygosity file
    het_df = pd.read_csv(f"{prefix}.het", sep='\t')
    het_df = het_df.rename(columns={'IID': 'Sample_ID'})

    # Combine files
    merged_df = miss_df.merge(het_df, on='Sample_ID', how='outer', suffixes=['_miss', '_het'])

    # Create heterozygosity columns 
    merged_df['heterozygosity'] = (merged_df['OBS_CT_het'] - merged_df['O(HOM)'])/merged_df['OBS_CT_het']
            
    # Calculate outlier limits for gbs samples
    gbs_het = merged_df.query('Sample_ID.str.contains("GBS")')['heterozygosity']
    gbs_het_lower_limit, gbs_het_upper_limit = outlier_limits(series=gbs_het, n=4)
    
    # Calculate outlier limits for riptide samples
    riptide_het = merged_df.query('~Sample_ID.str.contains("GBS")')['heterozygosity']
    riptide_het_lower_limit, riptide_het_upper_limit = outlier_limits(series=riptide_het, n=5)
    
    # Make outlier limits into columns
    merged_df['het_lower_limit'] = merged_df.apply(lambda row: gbs_het_lower_limit if 'GBS' in row['Sample_ID'] else riptide_het_lower_limit, axis=1)
    merged_df['het_upper_limit'] = merged_df.apply(lambda row: gbs_het_upper_limit if 'GBS' in row['Sample_ID'] else riptide_het_upper_limit, axis=1)	      
        
    # Create a heterogzygosity outlier status column
    merged_df['hetqc_outcome'] = merged_df.apply(lambda row: het_outlier(het_value=row['heterozygosity'],
                                                                        lower_limit=row['het_lower_limit'], 
                                                                        upper_limit=row['het_upper_limit'], 
                                                                        miss_value=row['missingness_perc'], 
                                                                        miss_limit=10), 
                                                                        axis=1)  
                                                                                    
    # Export dataframe
    merged_df[['Sample_ID', 'heterozygosity', 'het_lower_limit', 'het_upper_limit', 'hetqc_outcome']].to_csv(f"{prefix}.processed_het", index=False)


##############################################################################################################################################################################

def calculate_sample_pca_outliers(prefix, metadata_filepath, dtype_mapping):

    # Create a file containig outlier samples for exclusion

        # Read in pca file
    pca = pd.read_csv(f"{prefix}.eigenvec", sep="\t")
    pcs = pca[['PC1', 'PC2']]

        # Compute mean and covariance
    mean_vec = pcs.mean().values
    cov_matrix = np.cov(pcs.T)
    inv_cov_matrix = np.linalg.inv(cov_matrix)

        # Mahalanobis distance for each point
    pca['mahalanobis'] = pcs.apply(lambda x: mahalanobis(x, mean_vec, inv_cov_matrix), axis=1)

        # Define outliers: e.g., distance > 3 standard deviations
    threshold = pca['mahalanobis'].mean() + 3 * pca['mahalanobis'].std()

        # Create column indicating whether a SNP is an outlier
    pca['mahalanobis_outlier'] = pca['mahalanobis'].apply(lambda x: 'yes' if x > threshold else 'no')

        # Create a file of the outlier samples to exclude from PCA rerun
    outlier_samples = pca.query('mahalanobis_outlier == "yes"')['IID']
    outlier_samples.to_csv(f"{prefix}_outlier_samples.txt", header=['#IID'], index=False)
    print(f"Finished writing out outlier samples to {prefix}_outlier_samples.txt")
    
    # Create a file containing archival samples for exclusion

        # Read in pca file
    multi_flowcell_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping)

        # Query for archival samples
    archival_samples = multi_flowcell_df.query('Sample_Project.str.contains("archival|woods_rnaseq|retired_breeders")')['Sample_ID']

        # Combine archival and outlier samples
    archival_outlier_samples_set = set(outlier_samples.tolist() + archival_samples.tolist())

        # Create a dataframe from the archival and outlier samples for export
    archival_outlier_samples_df = pd.DataFrame({'$IID': list(archival_outlier_samples_set)})
    archival_outlier_samples_df.to_csv(f"{prefix}_archival_and_outlier_samples.txt", index=False)
    print(f"Finished writing out archival, outlier samples to {prefix}_archival_and_outlier_samples.txt")
