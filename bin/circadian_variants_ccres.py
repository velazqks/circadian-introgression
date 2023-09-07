#!/usr/bin/env python
# circadian_variants_ccres.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2022/03/25
#
# Description: Identify candidate cis-regulatory elements (cCREs) in the flanking regions 
#              near circadian genes, within 1Mb.
#              Source: https://doi.org/10.1038/s41586-020-2493-4
#
"""""""""""""""""""""""""""""""""


# DATA FILES
CCRES_FILE = '../data/raw_cCREs.liftOver.to.Hg19.bed.xz'
SNP_FILE = '../data/circadian_variants.bed.xz'

# OUTPUT FILE
OUTPUT_FILE = '../data/circadian_variants_ccres.bed'


import pybedtools
import pandas as pd


def intersect_a_and_b(df_a,df_b):
    df_a = df_a.add_suffix('_wa')
    df_b = df_b.add_suffix('_wb')
    df_a_cols = df_a.columns.values.tolist()
    df_b_cols = df_b.columns.values.tolist()
    a = pybedtools.BedTool.from_dataframe(df_a)
    b = pybedtools.BedTool.from_dataframe(df_b)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=df_a_cols+df_b_cols)
    return a_and_b_df


# LOAD DATA
ccres = pd.read_csv(CCRES_FILE, compression='xz', sep='\t')
snps = pd.read_csv(SNP_FILE, compression='xz', sep='\t')

# INTERSECT CIRCADIAN FLANKING SNPS AND CCRES
circadian_ccres = intersect_a_and_b(snps,ccres)

# FILTER COLUMNS
circadian_ccres = circadian_ccres.iloc[:,[0,1,2,3,4,5,6,7,-1]]

# REMOVE SUBSTRING FROM COLUMN NAMES
circadian_ccres.columns = circadian_ccres.columns.str.replace('_wa|_wb', '')


# SAVE DATA
circadian_ccres.to_csv(OUTPUT_FILE+'.xz', sep='\t', compression='xz', index=False)
