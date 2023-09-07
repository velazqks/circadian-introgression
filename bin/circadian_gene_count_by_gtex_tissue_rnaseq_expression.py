#!/usr/bin/env python
# circadian_gene_count_by_gtex_tissue_rnaseq_expression.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: To understand how the difference in power in each GTEx tissue affects 
#              the significance values of tissue-specific expression enrichment analysis, 
#              we test the fraction of circadian genes expressed in each GTEx tissue.
#              Expressed genes are defined as genes containing an RNASeq equal or above 1.
#
"""""""""""""""""""""""""""""""""


# INPUT DATA
CIRCADIAN_FILE = '../data/circadian_genes.list'
GTEX_RNASEQ_FILE = '../data/raw_GTEx_v8_RNASeQCv1.1.9_gene_median_tpm.gct.xz'

# OUTPUT DATA
OUTPUT_FILE = '../data/circadian_gene_count_by_gtex_tissue_rnaseq_expression.tsv'


import pandas as pd


# LOAD DATA
circadian_genes = pd.read_csv(CIRCADIAN_FILE, sep='\t')
rnaseq = pd.read_csv(GTEX_RNASEQ_FILE, sep='\t', compression='xz')
rnaseq.rename(columns={'Name':'GeneID'}, inplace=True)
# GTEx RNASeq has some outdated gene names. We will use GeneName from the circadian gene set
rnaseq.drop('Description', axis=1, inplace=True)


# FILTER THE GTEX RNASEQ DATAFRAME TO EXTRACT CIRCADIAN GENES
circadian_gtex_rnaseq = pd.merge(circadian_genes['GeneID'], rnaseq, on=['GeneID'])

# Set GeneID as index and then transpose the dataframe with GeneID as header
circadian_gtex_rnaseq_t = circadian_gtex_rnaseq.set_index(circadian_gtex_rnaseq.columns[0]).T

# RENAME TISSUES TO MATCH THE FORMAT IN THE GTEX V8 EQTL DATASET
circadian_gtex_rnaseq_t.index = circadian_gtex_rnaseq_t.index.str.replace(' - | ','_',regex=True).str.replace('\(|\)','',regex=True)

# Extract expression values higher than 1
gtex_expression = pd.DataFrame(circadian_gtex_rnaseq_t[(circadian_gtex_rnaseq_t>1)].count(axis=1).reset_index())
gtex_expression.columns = ['Tissue','Expressed_genes']
gtex_expression['Expressed_genes_percent'] = gtex_expression['Expressed_genes']/len(circadian_gtex_rnaseq)
gtex_expression = gtex_expression.sort_values(by='Expressed_genes', ascending=False).reset_index(drop=True)


# SAVE DATA
gtex_expression.to_csv(OUTPUT_FILE, sep='\t', index=True)
