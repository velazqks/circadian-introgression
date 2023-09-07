#!/usr/bin/env python
# circadian_predixcan_pctl.py
"""""""""""""""""""""""""""
# Keila Velazquez-Arcelay
#
# Description: Generate percentiles for where the archaic PrediXcan predictions fall 
#              in the human distribution
# 
"""""""""""""""""""""""""""


# INPUT FILES
CIRCADIAN_PRED = '../data/circadian_genes_predixcan_expression.tsv.xz'

# OUTPUT FILE
OUTPUT_FILE = '../data/circadian_genes_predixcan_pctl.tsv'


import pandas as pd
from scipy import stats


def generate_percentiles(df):
    # FIND THE PERCENTILE FOR THE ARCHAIC PREDICTIONS IN THE AMH DIST OF PREDICTIONS
    amh_ind = df.iloc[:,6:]
    l = []

    for i in range(len(df)):
        altai = stats.percentileofscore(amh_ind.iloc[i:i+1,:].values.tolist()[0],
                                        df['AltaiNeandertal'][i])
        vindija = stats.percentileofscore(amh_ind.iloc[i:i+1,:].values.tolist()[0],
                                          df['Vindija33.19'][i])
        denisova = stats.percentileofscore(amh_ind.iloc[i:i+1,:].values.tolist()[0],
                                          df['Denisova'][i])
        tissue = df['GTEx_Tissue'][i]
        geneid = df['GeneID'][i]
        genename = df['GeneName'][i]
        l.append(f'{geneid} {genename} {tissue} {altai} {vindija} {denisova}')
        
    # CREATE NEW DF
    df = pd.DataFrame([n.split() for n in l],
             columns=['GeneID','GeneName','GTEx_Tissue','AltaiPCTL','VindijaPCTL','DenisovaPCTL'])

    return df


# LOAD DATA
circadian_pred = pd.read_csv(CIRCADIAN_PRED, sep='\t', compression='xz')


# GENERATE PERCENTILES
circadian_pctl = generate_percentiles(circadian_pred)


# SAVE DATA
circadian_pctl.to_csv(OUTPUT_FILE,index=False,sep='\t')
