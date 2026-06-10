#!/usr/bin/env python
import pandas as pd
import os
import warnings
import sys
warnings.simplefilter(action='ignore', category=FutureWarning)
from utils.print_aesthetics import *

##############################################################################################################################################################################

def bad_snps(plink_prefix, novel_unknown_snps_filepath, chromvar):

    # Load in files
    # merdf = f"{plink_prefix}_mer_snp.csv"
    noveldf = pd.read_csv(f"{novel_unknown_snps_filepath}/chr{chromvar}_discordant_novel_snps_to_drop.csv", header=None, names=['chrom', 'pos', 'id'])
    infodf = pd.read_csv(f"{plink_prefix}.info", sep='\t').rename(columns={'[3]CHROM:[4]POS': 'id', '[5]INFO_SCORE': 'info'})
    if chromvar == 'X':
        infodf['info'] = infodf['info'].str.split(',').apply(lambda x: sum(float(item) for item in x) / 2)

    # Obtain list of SNPs to filter out
    bad_info = infodf.query('info < 0.9')['id'].tolist()
    bad_novel = noveldf['id'].tolist()
    bad_snps = list(set(bad_info + bad_novel))
    
    # Export bad SNPs list
    expdf = pd.DataFrame(bad_snps, columns=['id'])
    expdf[['chrom', 'pos']] = expdf['id'].str.split(':', expand=True)
    expdf['pos'] = expdf['pos'].astype(int)
    expdf = expdf.sort_values(by='pos').reset_index(drop=True)
    expdf.to_csv(f"{plink_prefix}_bad_snps.tsv", columns=['id'], sep='\t', header=False, index=None)

