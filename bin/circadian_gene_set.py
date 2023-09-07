#!/usr/bin/env python
# circadian_gene_set.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2023
#
# Description: Imports lists of circadian genes from different sources of evidence.
# Generates sets of circadian genes at different confidence levels based on the 
# intersection with the sources of evidence:
#     High: Genes expertly curated by Dr. McMahon or genes from 3/4 sources
#     Medium: Evidence from 2/4 sources
#     Low: Evidenced only in 1 source
# 
# https://pathcards.genecards.org/Card/circadian_clock?queryString=circadian
#
"""""""""""""""""""""""""""

# LARGE FILES
GENCODE_FILE = '../data/raw_gencode.v29lift37.bed'

# MY DATA
BIOSYS_FILE = '../data/raw_circadian_genes_biosystems.tsv'
CGDB_FILE = '../data/raw_circadian_genes_cgdb_experimental.tsv'
GO_FILE = '../data/raw_circadian_genes_go.tsv'
GWAS_FILE = '../data/raw_circadian_genes_gwas.tsv'
MCMAHON_FILE = '../data/raw_circadian_genes_mcmahon.tsv'

# OUTPUT FILE
OUTPUT = '../data/circadian_genes_candidate.tsv'


import pandas as pd
from functools import reduce


# Turn a list of genes from a source of evidence to a df
def list_to_column(l,source):
    l = l.rename(None)
    df = pd.DataFrame(l, columns=['GeneName'])
    df[source] = True
    df = pd.merge(gencode[['GeneID','GeneName']], df, on='GeneName')
    return df.drop_duplicates()


# Import Gencode df
try:
    gencode = pd.read_csv(GENCODE_FILE, sep='\t', encoding='utf-8')
except pd.errors.ParserError as e:
    print("ParserError:", e)
    
gencode.rename(columns={gencode.columns[3]:'GeneID',
                        gencode.columns[4]:'GeneName',
                        gencode.columns[5]:'GeneType'},
              inplace=True)

# Import circadian gene sources of evidence
biosystems = list_to_column(pd.read_csv(BIOSYS_FILE, sep='\t')['Symbol'],'Biosystems')
cgdb = list_to_column(pd.read_csv(CGDB_FILE, sep='\t')['GeneName'], 'CGDB')
go = list_to_column(pd.read_csv(GO_FILE, sep='\t')['GeneName'], 'GO')
gwas = list_to_column(pd.read_csv(GWAS_FILE, sep='\t')['GeneName'], 'GWAS')
mcmahon = list_to_column(pd.read_csv(MCMAHON_FILE, sep='\t')['GeneName'], 'McMahon')

# Merge sources of evidence in a membership table
sources = [biosystems,cgdb,go,gwas,mcmahon]
table = reduce(lambda left,right: pd.merge(left,right,on=['GeneID','GeneName'], 
                        how='outer'), sources).drop_duplicates()

# Insert empty Confidence column
table.insert(2, 'Confidence', None)

# Convert membership values from boolean to 1's and 0's
table.iloc[:,3:] = table.iloc[:,2:].fillna(0).astype(int)

# Create Confidence categories
table.loc[(table.iloc[:, 3:].sum(axis=1) == 3),'Confidence']='High'
table.loc[(table.iloc[:, 3:].sum(axis=1) == 2),'Confidence']='Medium'
table.loc[(table.iloc[:, 3:].sum(axis=1) == 1),'Confidence']='Low'
table.loc[(table['McMahon']==1),'Confidence']='High'

table.sort_values('GeneName', inplace=True)


# SAVE DATA
#table = table.fillna('NaN')
table.to_csv(OUTPUT, index=False, sep='\t')
