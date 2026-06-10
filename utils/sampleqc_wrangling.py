#!/usr/bin/env python
import pandas as pd
import os
import warnings
import sys
warnings.simplefilter(action='ignore', category=FutureWarning)
from utils.print_aesthetics import *

##############################################################################################################################################################################

def bad_samples(snp_filtered_outdir, current_round):

    chrom_list = list(range(1, 21)) + ['X', 'Y', 'M']   

    # Obtain missingness 
    missdf = {}

    for chromvar in chrom_list:
        tempdf = pd.read_csv(f"{snp_filtered_outdir}/chr{chromvar}_allsexes/{current_round}_chr{chromvar}_allsexes_snp_filtered.smiss", sep='\t')
        missdf[chromvar] = tempdf

    tempdf = pd.concat(missdf).reset_index(names=['chrom', 'old_index']).drop(columns=['old_index', 'F_MISS'])

    tempdf2 = tempdf.groupby('#IID')[['MISSING_CT', 'OBS_CT']].sum()
    tempdf2['F_MISS'] = tempdf2['MISSING_CT'] / tempdf2['OBS_CT']


    bad_samples = tempdf2.query('F_MISS > 0.1').reset_index()

    # Export bad SNPs list
    for chromvar in chrom_list:
        bad_samples.to_csv(f"{snp_filtered_outdir}/chr{chromvar}_allsexes/{current_round}_chr{chromvar}_allsexes_snp_filtered_bad_samples.tsv", columns=['#IID'], sep='\t', header=False, index=None)