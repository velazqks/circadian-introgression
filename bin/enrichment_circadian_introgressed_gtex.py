#!/usr/bin/env python
# enrichment_circadian_introgressed_gtex.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Hypergeometric distribution of introgressed variants in circadian genes 
#              with evidence of being eQTL in GTEx.
#
#              Introgressed variants are analyzed in two sets: Variants identified by Browning et al., 2018 
#              or variants identified by Browning and at least one other method.
# 
"""""""""""""""""""""""""""""""""


# INPUT DATA 
BG_LEN = 4272964    # Bi-allelic GTEx v8 eQTLs
CIRCADIAN_FILE = '../data/circadian_variants.bed.xz' # X  
CIRCADIAN_GTEX_FILE = '../data/circadian_variants_eqtl.bed.xz' # X  
BROWNING_FILE = '../data/raw_introgressed_browning2018.bed.xz'
BROWNING_GTEX_FILE = '../data/raw_GTEx_v8_browning2018.bed.xz'
INTROGRESSION_MAPS_FILE = '../data/introgression_maps.bed.xz'
# OUTPUT DATA 
OUTPUT_FILE = '../data/enrichment_circadian_introgressed_gtex.txt'


import os, sys
import pandas as pd
import numpy as np
import pybedtools
from scipy import stats


def load_file(file):
    df = pd.read_csv(file, sep='\t', compression='xz').iloc[:, :3].drop_duplicates()
    return df

def intersect_sets(df_a, df_b):
    df_a = df_a.add_suffix('_wa')
    df_b = df_b.add_suffix('_wb')
    a_cols = df_a.columns.values.tolist()
    b_cols = df_b.columns.values.tolist()
    a = pybedtools.BedTool.from_dataframe(df_a)
    b = pybedtools.BedTool.from_dataframe(df_b)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=a_cols + b_cols)
    return a_and_b_df

def fishers_exact(U,right,down,x,introgression_set):
    a = len(x)
    b = len(right) - a
    c = len(down) - a
    d = U - (a + b + c)

    #table = [[a, b], [c, d]]
    ar = np.array([[a, b],[c, d]])

    oddsratio, pvalue = stats.fisher_exact(ar)

    results = f'\nCircadian eQTLs in {introgression_set}\nOR: {oddsratio} \nP-VAL: {pvalue} \n{ar}'
        
    return results


# LOAD FILES
circadian = load_file(CIRCADIAN_FILE)
circadian_eqtl = load_file(CIRCADIAN_GTEX_FILE)
introgressed = load_file(BROWNING_FILE)
introgressed_eqtl = load_file(BROWNING_GTEX_FILE)

# FIND CIRCADIAN EQTLS
x_axis = pd.merge(circadian_eqtl,circadian)

# FIND INTROGRESSED EQTLS
y_axis = pd.merge(introgressed_eqtl,introgressed)

# IMPORT LIST OF INTROGRESSED SNPS DETECTED BY BROWNING2019 AND AT LEAST ONE OTHER METHOD
introgression_map = pd.read_csv(INTROGRESSION_MAPS_FILE, sep='\t', compression='xz')
introgression_map = introgression_map[introgression_map.iloc[:,4:].sum(axis=1) >= 2].iloc[:,:3].drop_duplicates()
y_axis_strict = pd.merge(introgression_map,introgressed_eqtl)

# FIND Q2: CIRCADIAN INTROGRESSED EQTLS
q2 = intersect_sets(x_axis,y_axis).drop_duplicates()
q2_strict = intersect_sets(x_axis,y_axis_strict).drop_duplicates()


# FISHER'S EXACT TEST
browning18 = fishers_exact(BG_LEN,x_axis,y_axis,q2,'Browning2018')
print(browning18)
browning18_plus = fishers_exact(BG_LEN,x_axis,y_axis_strict,q2,'Browning2018 and at least one other method')
print(browning18_plus)


# SAVE RESULTS
with open(OUTPUT_FILE, 'w') as f:
    f.write(browning18 + browning18_plus)
