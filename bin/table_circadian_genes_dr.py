#!/usr/bin/env python
# table_circadian_genes_dr.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Create a table containing a list of divergently regulated genes
#              in each of the archaic hominins. 
#              Input: Table 1.
#
"""""""""""""""""""""""""""""""""


# INPUT DATA
CIRCADIAN_FILE = '../data/table_circadian_genes.tsv'
ARCHAICS_FILE = '../data/circadian_genes_predixcan_dr.tsv'

# OUTPUT DATA
OUTPUT_FILE = '../data/table_circadian_genes_dr.xlsx'


import pandas as pd
#from functools import reduce


def load_archaics(file,ind):
    df = pd.read_csv(file, sep='\t')
    df = df[df[ind]==1]
    return df


# LOAD DATA
circadian_genes = pd.read_csv(CIRCADIAN_FILE, sep='\t').iloc[:,:3]
#circadian_genes = pd.read_excel('data/table_1.xlsx', engine='openpyxl').iloc[:,:3]
archaics = pd.read_csv(ARCHAICS_FILE, sep='\t')

# MERGE CIRCADIAN GENE DESCRIPTIONS AND ARCHAIC DIVERGENTLY REGULATED GENES
table = pd.merge(circadian_genes,archaics,on=['GeneID','GeneName']).sort_values(by='GeneName')

# SET MULTI-INDEX
table = table.set_index(['GeneID', 'GeneName', 'Description', 'GTEx_Tissue'])


# SAVE TO EXCEL
table.to_excel(OUTPUT_FILE, index=True)
