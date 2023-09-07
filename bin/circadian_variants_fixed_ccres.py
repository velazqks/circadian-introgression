#!/usr/bin/env python
# circadian_variants_fixed_ccres.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2022/06/09
#
# Description: Extract human- and archaic-specific circadian variants that are 
#              candidate cis-regulatory elements (cCREs).
#              Set of variants from Kuhlwilm and Boeckx 2019.
# 
"""""""""""""""""""""""""""""""""


# DATE FILES
HHMC_FILE = '../data/circadian_variants_hhmc.bed'
AHMC_FILE = '../data/circadian_variants_ahmc.bed'
CCRE_FILE = '../data/raw_cCREs.liftOver.to.Hg19.bed.xz'

# OUTPUT FILES
OUTPUT_AHMC = '../data/circadian_variants_ahmc_ccres.bed'
OUTPUT_HHMC = '../data/circadian_variants_hhmc_ccres.bed'

import pybedtools
import pandas as pd


def get_fixed_ccres(fixed,ccre):
    df = intersect_a_and_b(fixed,ccre)
    df = df.iloc[:,[0,1,2,3,4,5,6,7,-1]]
    df.columns = df.columns.str.replace('_wa|_wb', '')
    return df

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


# LOAD FILES
ahmc = pd.read_csv(AHMC_FILE, sep='\t')
hhmc = pd.read_csv(HHMC_FILE, sep='\t')
ccres = pd.read_csv(CCRE_FILE, sep='\t', compression='xz')

# FIND FIXED CIRCADIAN CCRES
ahmc_ccres = get_fixed_ccres(ahmc,ccres)
hhmc_ccres = get_fixed_ccres(hhmc,ccres)


# SAVE DATA
ahmc_ccres.to_csv(OUTPUT_AHMC, index=False, sep='\t')
hhmc_ccres.to_csv(OUTPUT_HHMC, index=False, sep='\t')
