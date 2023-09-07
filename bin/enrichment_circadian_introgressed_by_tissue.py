#!/usr/bin/env python
# enrichment_circadian_introgressed_by_tissue.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Find the probability of success of circadian eQTLs in each GTEx tissue,
#              in a background of all introgressed GTEx eQTLs.
#
"""""""""""""""""""""""""""""""""


# INPUT DATA
BROWNING_EQTL = '../data/raw_GTEx_v8_browning2018.bed.xz'
MAPS = '../data/introgression_maps.bed.xz'
C_GENES = '../data/circadian_genes.list'

# OUTPUT DATA
OUTPUT_FILE = '../data/enrichment_circadian_introgressed_by_tissue.tsv'


import pandas as pd
import numpy as np
import math
import pybedtools
from scipy import stats


def intersect_wa_wb(df_a, df_b):
    df_a = df_a.add_suffix('_wa')
    df_b = df_b.add_suffix('_wb')
    a_cols = df_a.columns.values.tolist()
    b_cols = df_b.columns.values.tolist()
    a = pybedtools.BedTool.from_dataframe(df_a)
    b = pybedtools.BedTool.from_dataframe(df_b)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=a_cols + b_cols)
    return a_and_b_df

def enrichment_in_tissues(U,right,down,x,t):
    # Fishers exact test
    a = len(x)
    b = len(right) - a
    c = len(down) - a
    d = (len(U) - (a + b + c))

    #table = [[a, b], [c, d]]
    ar = np.array([[a, b],[c, d]])

    oddsratio, pvalue = stats.fisher_exact(ar)
    
    return f'{t}, {oddsratio}, {pvalue}, {a}'


# LOAD FILES
c_genes = pd.read_csv(C_GENES,sep='\t').iloc[:,:1]
browning_eqtls = pd.read_csv(BROWNING_EQTL, sep='\t', compression='xz')
maps = pd.read_csv(MAPS, sep='\t', compression='xz')

# EXTRACT VARIANTS WITH EVIDENCE OF BEING INTROGRESSED IN BROWNING2018 AND AT LEAST ONE OTHER METHOD
filtered_maps = maps[(maps['Browning2018']==1) & (maps.iloc[:,5:].sum(axis=1)>=1)].iloc[:,:3].drop_duplicates()

# GET INTROGRESSED EQTLS
introgressed_snps = pd.merge(filtered_maps,browning_eqtls, on=['Chr','Start','End'])

# CREATE LIST OF TISSUES
tissues = introgressed_snps['Tissue'].drop_duplicates().values.tolist()

# DEFINE SAMPLES
# Y axis: Introgressed GTEx eQTLs matched to a tissue
y_axis = introgressed_snps.iloc[:,[0,1,2,-1]].drop_duplicates()
# Background: Unique introgressed eQTLs
background = introgressed_snps.iloc[:,0:3].drop_duplicates()
# X axis: Circadian SNPs that are introgressed eQTLs
x_axis = pd.merge(introgressed_snps,c_genes,left_on='GeneID',right_on='GeneID').iloc[:,[0,1,2]].drop_duplicates()

# PERFORM FISHER'S TEST ON EACH GTEX TISSUE
LIST = []
for i in tissues:
    # Select specific tissue
    tissue_iter = y_axis[y_axis['Tissue'].isin([i])].iloc[:, 0:3].drop_duplicates()
    # Find x: SNPs in tissue that are introgressed eQTLs
    x = intersect_wa_wb(x_axis,tissue_iter).iloc[:, 0:3].drop_duplicates()
    # Test for circadian introgressed enrichment in each GTEx tissue
    LIST.append( enrichment_in_tissues(background,x_axis,tissue_iter,x,i) )

# CREATE DATAFRAME FROM LIST
df = pd.DataFrame([n.split(', ') for n in LIST], columns=['Tissue','OR','P-value','Variants'])
df = df.astype({'Tissue': str, 'OR': float, 'P-value': float, 'Variants': str})

# REMOVE TISSUE WITH NO CIRCADIAN INTROGRESSED EQTLS
df = df[df['Tissue']!='Kidney_Cortex']

# CREATE OR LOG COLUMN
df['log2(odds ratio)'] = [math.log2(n) for n in df['OR']]

# ADD VARIANT COUNTS TO TISSUE NAMES
df['Tissues'] = df['Tissue'] + ' (' + df['Variants'] + ')'

# BONFERRONI CORRECTION
def myfunc(cat, pval):
    if pval>0.05:
        cat='>0.05'
    elif pval<=0.05 and pval>bonfe:
        cat='<=0.05'
    elif pval<=bonfe:
        cat='<='+str(round(bonfe,5))
    return cat

bonfe = 0.05/len(tissues)   # 49 tissues
df['Bonferroni'] = ''
df['Bonferroni'] = df.apply(lambda x: myfunc(x['Bonferroni'], x['P-value']), axis=1)


# SAVE DATA
df.to_csv(OUTPUT_FILE, sep='\t', index=False)
