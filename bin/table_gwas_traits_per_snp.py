#!/usr/bin/env python
# table_gwas_traits_per_snp.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Generate counts of GWAS phenotypes per introgressed SNP in two sets of SNPs:
#              Circadian and Non-circadian. The GWAS associations were extracted from
#              Opentargets.
#
"""""""""""""""""""""""""""""""""


# INPUT DATA
CIRCADIAN_INTROGRESSED = '../data/circadian_variants_introgressed.bed.xz'
CIRCADIAN_SNPS = '../data/circadian_variants.bed.xz'
INTROGRESSED = '../data/introgression_maps.bed.xz'
OPENTARGETS = '../data/opentargets_results_browning2018_signif_hg19.bed.xz'

# OUTPUT DATA
OUTPUT_FILE = '../data/plotting_gwas_traits_per_snp.tsv'


import pandas as pd


def traits_per_snp(assoc,pval,pheno):
    # Filter set of variants
    bonferroni = 0.00000005
    assoc = assoc[assoc[pval]<=bonferroni]
    assoc = assoc.drop(pval, axis=1).drop_duplicates()
    # Create locus column in the GWAS dataframe and merge with introgressed variants
    assoc['POS'] = assoc['Chr'] + '_' + assoc['End'].astype(str)
    index_i = pd.merge(assoc,introgressed,on='POS')['POS'].drop_duplicates()
    # Create two sets of variants: Circadian and non-circadian
    assoc_nc = assoc[~assoc['POS'].isin(circadian_i['POS'])]
    assoc_c = assoc[assoc['POS'].isin(circadian_i['POS'])]
    # Count phenotypes per variant
    b = assoc_nc.value_counts('POS').reset_index()#.reindex(index_i, fill_value=0).reset_index()
    a = assoc_c.value_counts('POS').reset_index()#.reindex(index_i, fill_value=0).reset_index()
    # Label the two sets
    b['SET'] = 'Non-Circadian'
    a['SET'] = 'Circadian'

    df = pd.concat([a,b])
    df.rename(columns={0:'COUNTS'}, inplace=True)
    
    return df


# LOAD DATA
# Opentargets associations
opentargets = pd.read_csv(OPENTARGETS, sep='\t', compression='xz')

# Introgressed variants
introgressed = pd.read_csv(INTROGRESSED, sep='\t', compression='xz')
# Keep only introgressed variants identified by Browning 2018 and at least one other method
introgressed = introgressed[introgressed.iloc[:,4:].sum(axis=1)>=2].iloc[:,:3].drop_duplicates()
# Create a Locus column
introgressed['POS'] = introgressed['Chr'].astype(str).str.cat(introgressed['End'].astype(str), sep='_')

# Circadian introgressed variants
circadian_i = pd.read_csv(CIRCADIAN_INTROGRESSED, 
                          sep='\t', compression='xz').iloc[:,:3].drop_duplicates()
circadian_i['POS'] = circadian_i['Chr'] + '_' + circadian_i['End'].astype(str)

# Extract introgressed variants in the Opentargets dataset
ot_introgressed = pd.merge(introgressed, opentargets, on=['Chr','Start','End'])

# Create counts of traits per SNP
ot_tps = traits_per_snp(ot_introgressed, 'pval', 'traitEfos')


ot_tps.to_csv(OUTPUT_FILE, sep='\t', index=False)