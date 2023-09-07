#!/usr/bin/env python
# enrichment_introgressed_gwas_or_no_gwas.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Identify if the sets of circadian and non-circadian introgressed variants 
#              have a significant difference in the proportion of variants with at 
#              least 1 GWAS association to variants with 0 associations.
# 
"""""""""""""""""""""""""""


# INPUT DATA
INTROGRESSED = '../data/raw_introgressed_browning2018.bed.xz'
OPENTARGETS = '../data/opentargets_results_browning2018_signif_hg19.bed.xz'
CIRCADIAN = '../data/circadian_variants.bed.xz'
THRESHOLD = 0.00000005

# OUTPUT DATA
OUTPUT_FILE = '../data/enrichment_introgressed_gwas_or_no_gwas.txt'


import pandas as pd
import numpy as np
import pybedtools
from scipy import stats


def load_data(filename):
    df = pd.read_csv(filename, sep='\t', compression='xz')
    return df


def create_locus_col(data):
    df = pd.DataFrame(data.iloc[:,0] + '_' + data.iloc[:,2].astype(str), 
                      columns=['Locus']).iloc[:,:3].drop_duplicates()
    return df


def fishers_exact(b,a,c,d):
    #table = [[a, b],[c, d]]
    ar = np.array([[len(a), len(b)],[len(c), len(d)]])
    
    oddsratio, pvalue = stats.fisher_exact(ar)
    
    results = f'OR: {oddsratio} \nP: {pvalue} \n{ar}\n'
    
    return results


# LOAD DATA
introgressed = load_data(INTROGRESSED)
introgressed_loci = create_locus_col(introgressed)

opentargets = load_data(OPENTARGETS)
opentargets = opentargets[opentargets['pval']<=THRESHOLD]
opentargets_loci = create_locus_col(opentargets)

circadian = load_data(CIRCADIAN)
circadian_loci = create_locus_col(circadian)


# DEFINE THE CIRCADIAN SET
introgressed_c = introgressed_loci[introgressed_loci['Locus'].isin(circadian_loci['Locus'])]
q2 = introgressed_c[introgressed_c['Locus'].isin(opentargets_loci['Locus'])]
q3 = introgressed_c[~introgressed_c['Locus'].isin(opentargets_loci['Locus'])]

# DEFINE THE NON-CIRCADIAN SET
introgressed_nc = introgressed_loci[~introgressed_loci['Locus'].isin(circadian_loci['Locus'])]
q1 = introgressed_nc[introgressed_nc['Locus'].isin(opentargets_loci['Locus'])]
q4 = introgressed_nc[~introgressed_nc['Locus'].isin(opentargets_loci['Locus'])]

# TEST FOR DIFFERENCE BETWEEN THE SETS
result = fishers_exact(q1,q2,q3,q4)
print(result)


# SAVE RESULTS
with open(OUTPUT_FILE, 'w') as f:
    f.write(result)
