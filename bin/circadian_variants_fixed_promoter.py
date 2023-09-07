#!/usr/bin/env python
# circadian_variants_fixed_promoter.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2022/06/09
#
# Description: Extract human- and archaic-specific variants in circadian promoter regions.
#              Variant set from Kuhlwilm and Boeckx 2019.
# 
"""""""""""""""""""""""""""""""""


# DATA FILES
PROMOTER_FILE = '../data/circadian_promoter.bed'
CIRCADIAN_GENES = '../data/circadian_genes.list'
AHMC_FILE = '../data/circadian_variants_ahmc.bed'
HHMC_FILE = '../data/circadian_variants_hhmc.bed'

# OUTPUT FILES
OUTPUT_AHMC = '../data/circadian_variants_ahmc_promoters.bed'
OUTPUT_HHMC = '../data/circadian_variants_hhmc_promoters.bed'


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

def find_circadian(df1,df2):
    df = intersect_a_and_b(df1,df2).iloc[:,:8]
    df.columns = df.columns.str.replace('_wa|_wb', '')
    return df


# LOAD CIRCADIAN PROMOTER REGIONS AND GENES
promoter = pd.read_csv(PROMOTER_FILE, sep='\t').iloc[:,:-1]

# LOAD SETS OF FIXED VARIANTS
circadian_ahmc = pd.read_csv(AHMC_FILE, sep='\t')
circadian_hhmc = pd.read_csv(HHMC_FILE, sep='\t')

# FIND FIXED CIRCADIAN VARIANTS IN PROMOTER REGIONS
promoter_ahmc = find_circadian(circadian_ahmc,promoter)
promoter_hhmc = find_circadian(circadian_hhmc,promoter)


# SAVE DATA
promoter_ahmc.to_csv(OUTPUT_AHMC, index=False, sep='\t')
promoter_hhmc.to_csv(OUTPUT_HHMC, index=False, sep='\t')
