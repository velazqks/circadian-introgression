#!/usr/bin/env python
# circadian_variants_fixed.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2022/06/09
#
# Description: Collects circadian variants from a set of human- and archaic-specific variants
#              published by Kuhlwilm and Boeckx 2019.
# 
"""""""""""""""""""""""""""""""""


# DATA FILES
HHMC_FILE = '../data/raw_kuhlwilm19_human_fixed_tony.bed'
AHMC_FILE = '../data/raw_kuhlwilm19_archaic_fixed.bed'
GENE_LOCI = '../data/circadian_genes.bed'
FLANKING_LOCI = '../data/circadian_genes_flanking.bed'

# OUTPUT DATA
OUTPUT_HHMC = '../data/circadian_variants_hhmc.bed'
OUTPUT_AHMC = '../data/circadian_variants_ahmc.bed'


import pandas as pd
import numpy as np
import pybedtools
#import warnings
#warnings.filterwarnings('ignore')


def load_kuhlwilm(filename):
    file = pd.read_csv(filename, sep='\t')
    file = file.iloc[:,[0,1,2,6,7]]
    return file

def load_circadian_gene_loci(filename, region):
    file = pd.read_csv(filename, sep='\t').iloc[:,:5]
    file['Region'] = region
    #file.rename(columns={'Start':'Start','End':'End'},inplace=True)
    return file

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

def wrangle_dfs(file):
    df = file.iloc[:,[0,1,2,3,4,-3,-2,-1]]
    df.columns = df.columns.str.replace(r'(_wa|_wb)', '')
    return df


# LOAD KUHLWILM VARIANTS
hhmc = load_kuhlwilm(HHMC_FILE)
ahmc = load_kuhlwilm(AHMC_FILE)

# LOAD CIRCADIAN LOCI
circadian_genes = load_circadian_gene_loci(GENE_LOCI, 'Gene')
circadian_flanking = load_circadian_gene_loci(FLANKING_LOCI, 'Flanking')

# CONCATENATE GENE AND FLANKING CIRCADIAN LOCI
circadian_loci = pd.concat([circadian_genes,circadian_flanking]).sort_values(by=['Chr','Start'])

# INTERSECT KUHLWILM VARIANTS AND CIRCADIAN LOCI
circadian_hhmc = intersect_a_and_b(hhmc,circadian_loci)
circadian_ahmc = intersect_a_and_b(ahmc,circadian_loci)

# SELECT RELEVANT COLUMNS AND REMOVE COLUMN NAME SUFFIXES
circadian_hhmc = wrangle_dfs(circadian_hhmc)
circadian_ahmc = wrangle_dfs(circadian_ahmc)


# SAVE DATA
circadian_hhmc.to_csv(OUTPUT_HHMC, index=False, sep='\t')
circadian_ahmc.to_csv(OUTPUT_AHMC, index=False, sep='\t')
