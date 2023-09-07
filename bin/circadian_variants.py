#!/usr/bin/env python
# circadian_variants.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2022/06/09
#
# Description: Filter the 1000 Genomes Project variants in circadian genes and flanking within 1Mb.
# 
"""""""""""""""""""""""""""""""""


# INPUT DATA FILES
GENOMESPRJ_FILE = '../general_data/KGP_phase3_v5c_20130502_sites.bed.xz'
GENES_FILE = '../data/circadian_genes.bed'
FLANKING_FILE = '../data/circadian_genes_flanking.bed'

# OUTPUT DATA FILES
OUTPUT_FILE = '../data/circadian_variants.bed'


import pandas as pd
import pybedtools


def load_circadian(file,region):
    df = pd.read_csv(file, sep='\t').iloc[:,:5]
    df['Region'] = region
    return df

def get_variants(x,y):
    x = x.add_suffix('_wa')
    y = y.add_suffix('_wb')
    x_cols = x.columns.values.tolist()
    y_cols = y.columns.values.tolist()
    
    a = pybedtools.BedTool.from_dataframe(x)
    b = pybedtools.BedTool.from_dataframe(y)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=x_cols+y_cols)
    
    return a_and_b_df


# LOAD DATA
genomesprj = pd.read_csv(GENOMESPRJ_FILE, sep='\t', compression='xz', low_memory=False)
circadian_genes = load_circadian(GENES_FILE,'Gene')
circadian_flanking = load_circadian(FLANKING_FILE,'Flanking')

# CONCATENATE ALL CIRCADIAN GENE LOCI
circadian_gene_loci = pd.concat([circadian_genes,circadian_flanking])
circadian_gene_loci.sort_values(by=['Chr','Start','GeneName'], inplace=True)

# MERGE CIRCADIAN GENE (AND FLANKING) LOCI AND 1KGP VARIANTS
# This line takes several minutes to process.
circadian_variants = get_variants(genomesprj,circadian_gene_loci)

# FILTER COLUMNS
circadian_variants = circadian_variants.iloc[:,[0,1,2,3,4,-3,-2,-1]].drop_duplicates()

# REMOVE COLUMN NAME SUFFIXES
circadian_variants.columns = circadian_variants.columns.str.replace('(_wa|_wb)', '')


# SAVE DATA IN lzma-COMPRESSED FORMAT
circadian_variants.to_csv(OUTPUT_FILE+'.xz', sep='\t', compression='xz', index=False)
