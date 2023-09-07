#!/usr/bin/env python
# table_circadian_genes_sav.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Create a table containing a list of archaic variants predicted to be  
#              splice-altering by SpliceAI. Webscrape the descriptions for the genes 
#              containing the SAVs. Input: circadian_variants_sav.tsv.
#
"""""""""""""""""""""""""""""""""


# INPUT DATA
CIRCADIAN_FILE = '../data/table_circadian_genes.tsv'
SAVS_FILE = '../data/circadian_variants_sav.tsv'

# OUTPUT DATA
OUTPUT_FILE = '../data/table_circadian_genes_sav.xlsx'


import pandas as pd


# LOAD DATA
circadian_genes = pd.read_csv(CIRCADIAN_FILE, sep='\t').iloc[:,:3]
savs = pd.read_csv(SAVS_FILE, sep='\t')

# ADD LOCUS COLUMN FOR EACH VARIANT
savs.insert(6,'Locus', savs['Chr'] + '_' + savs['End'].astype(str))
savs = savs.iloc[:,4:]

# MERGE GENE DESCRIPTION WITH VARIANT LOCI AND SAV INFORMATION FOR EACH ARCHAIC
table = pd.merge(circadian_genes, savs, on=['GeneID','GeneName']).sort_values(by='GeneName')

# SET MULTI-INDEX
table = table.set_index(['GeneID', 'GeneName', 'Description', 'Locus'])


# SAVE TO EXCEL
table.to_excel(OUTPUT_FILE, index=True)
