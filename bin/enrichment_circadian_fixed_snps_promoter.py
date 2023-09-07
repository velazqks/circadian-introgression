#!/usr/bin/env python
# enrichment_circadian_fixed_snps_promoter.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Enrichment analysis of human and archaic lineage-specific SNPs in circadian promoters.
#              Set of SNPs from Kuhlwilm and Boeckx, 2019. 
# 
"""""""""""""""""""""""""""""""""


# INPUT DATA 
KUHLWILM_FILE = '../data/raw_kuhlwilm2019_filtered.bed.xz'
HHMC = '../data/raw_kuhlwilm19_human_fixed_tony.bed'
AHMC = '../data/raw_kuhlwilm19_archaic_fixed.bed'
PROMOTER_FILE = '../data/circadian_promoter.bed'

# OUTPUT DATA
OUTPUT_FILE = '../data/enrichment_circadian_fixed_snps_promoter.txt'


import pandas as pd
import numpy as np
import pybedtools
from scipy import stats


def get_variant_set(filename):
    file = pd.read_csv(filename, sep='\t')
    universe = fix_pos_column(file)
    U = len(file['POS'].drop_duplicates())
    return universe, U

def fix_pos_column(data):
    df = pd.DataFrame(columns=['Chr','End'])
    df[['Chr','End']] = data['POS'].str.split(':', expand=True)
    df['Start'] = df['End'].map(int)-1
    df['End'] = df['End'].map(int)
    df['Chr'] = df['Chr'].replace('^','chr',regex=True)
    df = df[['Chr', 'Start', 'End']]
    return df

def get_x_axis(universe):
    # Import loci
    x = pd.read_csv(PROMOTER_FILE, sep='\t')
    # Intersect loci with set of SNPs
    snps = intersect_a_and_b(universe,x).iloc[:,:3].drop_duplicates()
    snps.columns = snps.columns.str.replace('_wa|_wb', '')
    snps_len = len(snps[['Chr','Start','End']].drop_duplicates())
    return snps, snps_len

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

def get_y_axis(filename):
    file = pd.read_csv(filename, sep='\t')
    y = fix_pos_column(file)
    return y
    
def fishers_exact(x,row_tot,col_tot,tot,group):
    a = len(x)
    b = len(row_tot) - a
    c = len(col_tot) - a
    d = tot - (a + b + c)

    #table = [[a, b],[c, d]]
    ar = np.array([[a, b],[c, d]])
    
    oddsratio, pvalue = stats.fisher_exact(ar)
    
    results = f'\n{group} specific enrichment \nOR: {oddsratio} \nP: {pvalue} \n{ar}\n'
    
    return results

def contingency_table(ar,group):
    df = pd.DataFrame(ar, columns=[f"{group} specific variants", f"Non-{group} specific variants"])
    df.index = ["Circadian variants", "Non-Circadian variants"] 
    df.loc['Column_Total']= df.sum(numeric_only=True, axis=0)
    df.loc[:,'Row_Total'] = df.sum(numeric_only=True, axis=1)
    return df


# LOAD GENERAL SET OF VARIANTS
all_snps, U = get_variant_set(KUHLWILM_FILE)

# FIND ALL CIRCADIAN SNPS IN THE SET OF SNPS
circadian_snps, circadian_len = get_x_axis(all_snps)

# LOAD FIXED VARIANTS
y_a = get_y_axis(AHMC)
y_h = get_y_axis(HHMC)

# FIND Q2: Fixed circadian
# Archaic
y_a_circadian = pd.merge(circadian_snps,y_a, on=['Chr','Start','End'])
# Human
y_h_circadian = pd.merge(circadian_snps,y_h, on=['Chr','Start','End'])

# FISHER'S EXACT TEST
# Archaic specific
aar = fishers_exact(y_a_circadian, circadian_snps, y_a, U, 'Archaic')
print(aar)
# Human specific
har = fishers_exact(y_h_circadian, circadian_snps, y_h, U, 'Human')
print(har)
# contingency_table(aar,'Archaic')
# contingency_table(har,'Human')


# SAVE RESULTS
with open(OUTPUT_FILE, 'w') as f:
    f.write(aar + har)
