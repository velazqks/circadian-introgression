#!/usr/bin/env python
# circadian_genes_divergently_regulated.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Find divergently regulated circadian genes, defined as PrediXcan
#              predictions with an empirical p-value = 0.
# 
#              PrediXcan results published in Colbran et al., 2019
#              https://github.com/colbrall/neanderthal_predixcan_manuscript
#
"""""""""""""""""""""""""""


# INPUT DATA
ALTAI_FILE = '../data/raw_predixcan_pvals_circadian_Altai.tsv'
VINDIJA_FILE = '../data/raw_predixcan_pvals_circadian_Vindija.tsv'
DENISOVA_FILE = '../data/raw_predixcan_pvals_circadian_Denisova.tsv'

# OUTPUT DATA
OUTPUT_FILE = '../data/circadian_genes_predixcan_dr.tsv'


import pandas as pd
from functools import reduce


def load_and_wrangle_predictions(file, ind):
    # Load archaic p-values file
    df = pd.read_csv(file, sep='\t')
    # Convert the tissue columns to rows using melt
    df = df.melt(id_vars=['GeneID','GeneName'], var_name='GTEx_Tissue', value_name='p-value')
    # Extract the significant predictions: p-value=0
    df = df[df['p-value']==0.0]
    # Keep only GeneID and GTEx_Tissue columns
    df = df[['GeneID', 'GeneName', 'GTEx_Tissue']]
    df = df.sort_values(by='GeneName')
    # Add a column for each archaic to indicate where it has significant predictions for each gene
    df[ind] = 1
    return df


# FIND SIGNIFICANTLY DIVERGENTLY REGULATED GENES
altai_signif = load_and_wrangle_predictions(ALTAI_FILE,'Altai')
vindija_signif = load_and_wrangle_predictions(VINDIJA_FILE,'Vindija')
denisova_signif = load_and_wrangle_predictions(DENISOVA_FILE,'Denisova')

# MERGE ALL ARCHAIC DFS
DFS = [altai_signif, vindija_signif, denisova_signif]
COLS = ['GeneID','GeneName','GTEx_Tissue']
archaics_signif = reduce(lambda left,right: pd.merge(left,right,on=COLS,how='outer'), DFS)
archaics_signif.iloc[:,-3:] = archaics_signif.iloc[:,-3:].fillna('0').astype(int)

# DR in all archaics
#print(archaics_signif[archaics_signif.sum(axis=1) == 3])


# SAVE FILE
archaics_signif.to_csv(OUTPUT_FILE, sep='\t', index=False)
