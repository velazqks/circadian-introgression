#!/usr/bin/env python
# circadian_gene_count_by_predixcan_tissue_prediction.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: To understand how the difference in power in each GTEx tissue affects 
#              the significance values of tissue-specific expression enrichment analysis, 
#              we test how many circadian genes contain PrediXcan predictions in each tissue.
#
"""""""""""""""""""""""""""""""""


# INPUT DATA
CIRCADIAN_FILE = '../data/circadian_genes.list'
PREDIXCAN_FILE = '../data/circadian_genes_predixcan_expression.tsv.xz'

# OUTPUT DATA
OUTPUT_FILE = '../data/circadian_gene_count_by_predixcan_tissue_prediction.tsv'


import pandas as pd


# LOAD DATA
circadian_genes = pd.read_csv(CIRCADIAN_FILE, sep='\t')
circadian_pred = pd.read_csv(PREDIXCAN_FILE, sep='\t', compression='xz').iloc[:,:3].drop_duplicates()


# COUNT NUMBER OF CIRCADIAN GENES WITH PREDIXCAN PREDICTIONS IN EACH TISSUE
predixcan = circadian_pred.value_counts('GTEx_Tissue').reset_index()
predixcan.columns = ['Tissue', 'Gene_Count']
predixcan['Gene_Count_Percentage'] = predixcan['Gene_Count'] / len(circadian_genes)


# SAVE DATA
predixcan.to_csv(OUTPUT_FILE, sep='\t', index=False)
