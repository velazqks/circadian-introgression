#!/usr/bin/env python
# enrichment_circadian_fixed_snps_regulatory.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Enrichment analysis of human and archaic lineage-specific SNPs in circadian 
#              regulatory regions (cCREs). Set of SNPs from Kuhlwilm and Boeckx, 2019.
# 
"""""""""""""""""""""""""""""""""

# INPUT DATA
KUHLWILM_FILE = '../data/raw_kuhlwilm2019_filtered.bed.xz'
HHMC = '../data/raw_kuhlwilm19_human_fixed_tony.bed'
AHMC = '../data/raw_kuhlwilm19_archaic_fixed.bed'
CCRE_FILE = '../data/raw_cCREs.liftOver.to.Hg19.bed.xz'
CIRCADIAN_FLANKING = '../data/circadian_genes_flanking.bed'

# OUTPUT DATA
OUTPUT_FILE = '../data/enrichment_circadian_fixed_snps_regulatory.txt'


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
    x = pd.read_csv(CIRCADIAN_FLANKING, sep='\t')
    # Import cCREs and find cCREs in the general set of SNPs
    ccres = load_ccres(universe)
    circadian_ccres = intersect_a_and_b(ccres,x).iloc[:,:3].drop_duplicates()
    circadian_ccres.columns = circadian_ccres.columns.str.replace('_wa|_wb', '')
    circadian_len = len(circadian_ccres[['Chr','Start','End']].drop_duplicates())
    return circadian_ccres

def load_ccres(u):
    ccres = pd.read_csv(CCRE_FILE, sep='\t')
    ccres_kuhlwilm = intersect_a_and_b(u,ccres).iloc[:,:3].drop_duplicates()
    ccres_kuhlwilm.columns = ccres_kuhlwilm.columns.str.replace('_wa|_wb', '')
    return ccres_kuhlwilm

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
    file = pd.read_csv(filename, sep='\t', low_memory=False)
    y = fix_pos_column(file)
    y_len = len(file['POS'].drop_duplicates())
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

# FIND ALL CIRCADIAN cCRES IN THE SET OF SNPS
circadian_ccres = get_x_axis(all_snps)

# LOAD FIXED VARIANTS
y_a = get_y_axis(AHMC)
y_h = get_y_axis(HHMC)

# FIND Q2: Fixed circadian SNPs
# Archaic
y_a_circadian = pd.merge(circadian_ccres, y_a, on=['Chr','Start','End'])
# Human
y_h_circadian = pd.merge(circadian_ccres, y_h, on=['Chr','Start','End'])

# FISHER'S: Archaic specific.
aar = fishers_exact(y_a_circadian, circadian_ccres, y_a, U, 'Archaic')
print(aar)
#print(contingency_table(aar,'Archaic'))

# FISHER'S: Human specific
har = fishers_exact(y_h_circadian, circadian_ccres, y_h, U, 'Human')
print(har)
#print(contingency_table(har,'Human'))


# SAVE RESULTS
with open(OUTPUT_FILE, 'w') as f:
    f.write(aar + har)
