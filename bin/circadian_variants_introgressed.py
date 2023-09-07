#!/usr/bin/env python
# circadian_variants_introgressed.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2022/06/09
#
# Description: Extract circadian variants from a set of introgressed variants
#              published by Browning et al., 2018. The set of circadian variants
#              includes variants inside circadian genes, cCREs flanking circadian 
#              genes within 1Mb, and promoter variants flanking the TSS by
#              -5kb upstream and 1kb downstream.
# 
"""""""""""""""""""""""""""""""""


# INPUT DATA
BROWNING = '../data/raw_introgressed_browning2018.bed.xz'
CIRCADIAN_SNPS = '../data/circadian_variants.bed.xz'
CCRE_SNPS = '../data/circadian_variants_ccres.bed.xz'
PROMOTER_SNPS = '../data/circadian_variants_promoter.bed.xz'

# OUTPUT DATA
OUTPUT_FILE = '../data/circadian_variants_introgressed.bed.xz'


import pandas as pd
import pybedtools

    
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

def collapse_columns(df,col):
    cols = df.columns.values.tolist()
    cols.remove(col)
    df[col] = df[col].astype(str)
    df = df.groupby(cols)[col].apply('|'.join).reset_index()
    return df


# LOAD FILES
browning = pd.read_csv(BROWNING, sep='\t', compression='xz').iloc[:,:5]
circadian_snps = pd.read_csv(CIRCADIAN_SNPS, sep='\t', compression='xz')#.iloc[:,[0,1,2,3,4,6]]
circadian_ccres = pd.read_csv(CCRE_SNPS, sep='\t', compression='xz').iloc[:,:-1]
circadian_promoters = pd.read_csv(PROMOTER_SNPS, sep='\t')

# 
circadian_genic = circadian_snps[circadian_snps['Region']=='Gene']
circadian_promoters['Region'] = 'Promoter'
circadian_ccres['Region'] = 'Regulatory'

# MERGE ALL CIRCADIAN SNP DATAFRAMES
dfs = pd.concat([circadian_genic,circadian_promoters,circadian_ccres])

# Collapse duplicated loci in multiple Regions
# dfs = collapse_columns(dfs,'Region')

# INTERSECT INTROGRESSED AND CIRCADIAN VARIANTS
circadian_introgressed = intersect_a_and_b(browning,dfs)

# FILTER COLUMNS
circadian_introgressed = circadian_introgressed.iloc[:,[0,1,2,3,4,-3,-2,-1]].drop_duplicates()

# REMOVE COLUMN SUFFIX
circadian_introgressed.columns = circadian_introgressed.columns.str.replace('_wa|_wb', '')


# SAVE DATA
circadian_introgressed.to_csv(OUTPUT_FILE, index=False, sep='\t', compression='xz')
