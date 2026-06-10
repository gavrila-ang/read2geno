#!/usr/bin/env python
import pandas as pd
import subprocess
import os
import sys
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=pd.errors.SettingWithCopyWarning)
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline



##############################################################################################################################################################################

def create_svm_input(metadata_filepath, samtools_view_filepath, exclude_sampleids_list, stage, seq_method, dtype_mapping, output_directory):
    
    # Read in the flowcell data
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    if seq_method == 'riptide':
        metadata_df = metadata_df.query(f"seq_method.isin(['riptide', 'LSW_lcwgs'])").reset_index(drop=True)
    elif seq_method == 'gbs':
        metadata_df = metadata_df.query(f"seq_method == 'gbs'").reset_index(drop=True)

    else:
        print(f"ERROR: Invalid sequencing method: {seq_method}")
        sys.exit(0)

    # Read in the samtools view concatenated file
    smtl_df = pd.read_csv(samtools_view_filepath)

    # Join metadata and samtools view concatenated file
    svm_data_full = smtl_df.merge(metadata_df[["Sample_ID", "sex", "Library_ID", "Sample_Project", "seq_method"]], how="inner", on="Sample_ID")

    # Remove samples from the exclude list
    svm_data_full = svm_data_full[~svm_data_full.Sample_ID.isin(exclude_sampleids_list)].reset_index(drop=True)

    # Set missing sex data as a category
    svm_data_full['sex'] = svm_data_full['sex'].fillna('missing')

    # Create normalised data
    svm_data_full['chrX'] = svm_data_full['chrX'].astype(int)
    svm_data_full['chrY'] = svm_data_full['chrY'].astype(int)
    svm_data_full['allchrom'] = svm_data_full['allchrom'].astype(int)
    svm_data_full['chrX_perc'] = ((svm_data_full.chrX / svm_data_full.allchrom) * 100).astype(float)
    svm_data_full['chrY_perc'] = ((svm_data_full.chrY / svm_data_full.allchrom) * 100).astype(float)

    # Drop all the samples with no sex recorded from the training data
    svm_data_dropna = svm_data_full[svm_data_full.sex!='missing'].reset_index(drop=True)
    
    # Export dataframe
    if stage == 'clean':
        svm_data_dropna.to_csv(f"{output_directory}/{stage}_{seq_method}_svm_train_data.csv", index=False)   
        svm_data_dropna.to_csv(f"{output_directory}/{stage}_{seq_method}_svm_predict_data.csv", index=False)                

    elif stage == 'impute':
        svm_data_full.to_csv(f"{output_directory}/{stage}_{seq_method}_svm_predict_data.csv", index=False)   
    
    else:
        print("ERROR: Invalid stage value.")
        sys.exit(1)


##############################################################################################################################################################################

def run_svm_linear(input_directory, dtype_mapping, stage, seq_method, output_directory, downsample='no'):
    
    # Read in the training data
    train_data = pd.read_csv(f"{input_directory}/{stage}_{seq_method}_svm_train_data.csv", dtype=dtype_mapping).reset_index(drop=True)
    
    # Read in the predict data
    predict_data = pd.read_csv(f"{input_directory}/{stage}_{seq_method}_svm_predict_data.csv", dtype=dtype_mapping).reset_index(drop=True)
    
    # Downsample training data
    if downsample != 'no':
        train_data = train_data.sample(frac=downsample)
    
    # Fit a linear SVM
    pipe = Pipeline([
        ("scale", StandardScaler()),
        ("model", SVC(kernel='linear', probability=True, random_state=98, max_iter=-1)
        )
    ])
    
    # Train model
    X_train = train_data[['chrX_perc', 'chrY_perc']]
    y_train = train_data['sex']
    model = pipe.fit(X_train, y_train)
    
    # Predict on model
    X_predict = predict_data[['chrX_perc', 'chrY_perc']]
    predicted_sex = model.predict(X_predict)

    # Turn predicted values into a dataframe
    predicted_sex_df = pd.DataFrame(predicted_sex, columns=['predicted_sex'])

    # Turn classification probabilities into a dataframe
    sex_proba = pipe.predict_proba(X_predict)
    sex_proba_df = pd.DataFrame(sex_proba, columns=['F_proba', 'M_proba'])

    # Concatenate the predictions, probabilities into the `predict data` dataset
    svm_pred_proba_df = pd.concat([predict_data, predicted_sex_df, sex_proba_df],
                                  axis=1)

    # Export dataframe
    svm_pred_proba_df.to_csv(f"{output_directory}/{stage}_{seq_method}_svm_output.csv", index=False)
    

##############################################################################################################################################################################

def create_sexqc_outcome_col(input_directory, dtype_mapping, stage, seq_method, gbs_sexqc_limit, riptide_sexqc_limit, output_directory):
    
    # Read in the svm data
    svm_data = pd.read_csv(f"{input_directory}/{stage}_{seq_method}_svm_output.csv", dtype=dtype_mapping).reset_index(drop=True)
    
    svm_data['sexqc_outcome'] = None

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

    for index, row in svm_data.iterrows():
        if (row['F_proba'] < undetermined_cutoff) and (row['M_proba'] < undetermined_cutoff):         
            svm_data.loc[index, 'sexqc_outcome'] = 'undetermined'
        elif row['sex'] == 'missing':
            svm_data.loc[index, 'sexqc_outcome'] = 'imputed'
        elif row['sex'] != row['predicted_sex']:
            svm_data.loc[index, 'sexqc_outcome'] = 'record_mismatch'
        elif row['sex'] == row['predicted_sex']:
            svm_data.loc[index, 'sexqc_outcome'] = 'kept_original'
        else:
            print(f"ERROR in assigning sexqc_outcome to {row.Sample_ID}.")
    
    # Export dataframe
    svm_data.to_csv(f"{output_directory}/{stage}_{seq_method}_svm_output_outcome.csv", index=False)


##############################################################################################################################################################################

def create_svm_clean_input(input_directory, dtype_mapping, stage, seq_method, output_directory):
    
    # Read in svm data
    svm_data = pd.read_csv(f"{input_directory}/{stage}_{seq_method}_svm_output_outcome.csv", dtype=dtype_mapping).reset_index(drop=True)

    if (stage == 'clean') and (seq_method == 'riptide'):
        clean_data = svm_data[(svm_data['sex']==svm_data['predicted_sex']) & ((svm_data['F_proba'] >= 0.8) | (svm_data['M_proba'] >= 0.8))]         
        clean_data.to_csv(f"{output_directory}/impute_{seq_method}_svm_train_data.csv", index=False)

    elif (stage == 'impute') and (seq_method == 'riptide'):
        print("SVM training data already cleaned.")

    elif (stage == 'clean') and (seq_method == 'gbs'):
        clean_data = svm_data[(svm_data['sex']==svm_data['predicted_sex']) & ((svm_data['F_proba'] >= 0.98) | (svm_data['M_proba'] >= 0.98))]         
        clean_data.to_csv(f"{output_directory}/impute_{seq_method}_svm_train_data.csv", index=False)

    elif (stage == 'impute') and (seq_method == 'gbs'):
        print("SVM training data already cleaned.")

    else:
        print("ERROR: Invalid stage and seq_method combination.")
        sys.exit(1)


##############################################################################################################################################################################

def create_plink_sex_data(stitch_samples_filepath, output_directory, current_round, dtype_mapping):

    # Read in the flowcell data    
    stitch_samples_df = pd.read_csv(stitch_samples_filepath, dtype=dtype_mapping).reset_index(drop=True)

    allsexes_sex_data = stitch_samples_df[['Sample_ID', 'predicted_sex']].rename(columns={'Sample_ID': '#IID', 'predicted_sex':'SEX'})
    allsexes_sex_data.to_csv(f"{output_directory}/{current_round}_plink_allsexes_sex_data.csv", header=True, sep='\t', index=False)
    
    female_sex_data = stitch_samples_df.loc[stitch_samples_df['predicted_sex'] == 'F', ['Sample_ID', 'predicted_sex']].rename(columns={'Sample_ID': '#IID', 'predicted_sex':'SEX'})
    female_sex_data.to_csv(f"{output_directory}/{current_round}_plink_female_sex_data.csv", header=True, sep='\t', index=False)
    
    male_sex_data = stitch_samples_df.loc[stitch_samples_df['predicted_sex'] == 'M', ['Sample_ID', 'predicted_sex']].rename(columns={'Sample_ID': '#IID', 'predicted_sex':'SEX'})
    male_sex_data.to_csv(f"{output_directory}/{current_round}_plink_male_sex_data.csv", header=True, sep='\t', index=False)    


##############################################################################################################################################################################

def svm_final_outputs(input_directory, metadata_filepath, dtype_mapping, current_round, skip_sexqc, output_directory):
    
    # Read in multi-flowcell sheet
    metadata_df = pd.read_csv(metadata_filepath, dtype=dtype_mapping).reset_index(drop=True)

    # Declare svm files
    gbs_file = f"{input_directory}/impute_gbs_svm_output_outcome.csv"
    riptide_file = f"{input_directory}/impute_riptide_svm_output_outcome.csv"
    
    # Concatenate svm outputs
    if os.path.exists(gbs_file):
        gbs_svm = pd.read_csv(gbs_file)
        riptide_svm = pd.read_csv(riptide_file)
        full_svm = pd.concat([gbs_svm, riptide_svm])
    else:
        full_svm = pd.read_csv(riptide_file)

    # Create STITCH samples
    if not skip_sexqc:
        stitch_data = full_svm[full_svm.sexqc_outcome.isin(['kept_original', 'imputed'])]
        stitch_data.to_csv(f"{output_directory}/{current_round}_stitch_samples.csv", index=False)
    else:
        full_svm.to_csv(f"{output_directory}/{current_round}_stitch_samples.csv", index=False)
    
    # Create a manual review file
    lab_check_data = full_svm.loc[full_svm.sexqc_outcome.isin(['undetermined', 'record_mismatch']),
                                  ['Sample_ID', 'Library_ID', 'Sample_Project', 'sex', 'predicted_sex', 
                                   'chrX_perc', 'chrY_perc', 'F_proba', 'M_proba', 'sexqc_outcome']]
    lab_check_data.to_csv(f"{output_directory}/manual_review_excluded_samples.csv", index=False)

    # Create add to exclude Sample ID
    full_svm = full_svm.merge(metadata_df[['Sample_ID', 'rfid', 'Sample_Barcode']], on='Sample_ID', how='left')
    lab_check_data = full_svm.query('sexqc_outcome.isin(["undetermined", "record_mismatch"])')
    lab_check_data['round_dropped'] = current_round
    lab_check_data['reason'] = 'step5_alignqc_sexqc' + '_' + lab_check_data['sexqc_outcome']
    lab_check_data['notes'] = None

    lab_check_data = lab_check_data[['Sample_ID', 'Library_ID', 'rfid', 'Sample_Barcode', 'round_dropped', 'reason', 'notes']]
    lab_check_data.to_csv(f"{output_directory}/add_to_exclude_sample_list.csv", index=False)      