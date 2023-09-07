#!/usr/bin/env python
# circadian_gene_loci.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2022/06/09
#
# Description: Create bed file of circadian genes and their 1Mb flanking loci.
#
"""""""""""""""""""""""""""""""""

# DATA FILES
GENCODE_FILE = '../data/raw_gencode.v29lift37.bed'
CIRCADIAN_FILE = '../data/circadian_genes_candidate.tsv'

# OUTPUT FILES
OUTPUT_GENES = '../data/circadian_genes.bed'
OUTPUT_FLANKING = '../data/circadian_genes_flanking.bed'
OUTPUT_GENE_LIST = '../data/circadian_genes.list'


import pandas as pd


def split_gene_columns(df,col):
    # Parse gencode file and reformat last column into different columns
    df = column_name_mapping(df)
    
    temp_df = df[8].str.split("; ", expand = True)[[0,col]] # , n = 2
    df["GeneID"] = temp_df[0]
    df["GeneName"] = temp_df[col]
    df.drop(columns=[8],inplace = True)
    
    df = replace_gene_strings(df)
    return df
    
def add_flanking_region(data):
    # Create flanking upstream region
    df_up = data.copy()
    df_up['End'] == df_up['Start']
    df_up['Start'] -= 1000000
    
    # Change value to 0 if Upstream limit is a negative value
    mask = df_up['Start'] < 0
    df_up.loc[mask,'Start'] = 0

    # Create a DataFrame with added values
    df_down = data.copy()
    df_down['Start'] == df_down['End']
    df_down['End'] += 1000000

    # Concatenate the original DataFrame with the subtracted and added DataFrames
    df = pd.concat([df_up, df_down], ignore_index=True)
    df.sort_values(by=['Chr', 'GeneName', 'Start'], inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


# LOAD GENCODE DATA
gencode = pd.read_csv(GENCODE_FILE, sep='\t')
gencode.rename(columns={gencode.columns[3]:'GeneID',
                        gencode.columns[4]:'GeneName',
                        gencode.columns[5]:'GeneType'},
              inplace=True)

# LOAD CIRCADIAN GENE SET
circadian = pd.read_csv(CIRCADIAN_FILE, sep='\t').iloc[:,:3]
circadian = circadian[~circadian['Confidence'].isin(['Low'])]

# ADD GENE LOCI
genic_region = pd.merge(gencode,circadian,on=['GeneID','GeneName'])
genic_region = genic_region[genic_region['Chr']!='chrY']

# ADD UPSTREAM AND DOWNSTREAM FLANKING LOCI
flanking_region = add_flanking_region(genic_region)

# EXTRACT LIST OF GENES
gene_list = genic_region[['GeneID','GeneName']].drop_duplicates().sort_values(by='GeneName')
# (head -n 1 circadian_genes.bed && tail -n +2 circadian_genes.bed | sort -k5) | cut -f4,5 >  circadian_genes.list


# SAVE DATA
genic_region.to_csv(OUTPUT_GENES, index=False, sep='\t')
flanking_region.to_csv(OUTPUT_FLANKING, index=False, sep='\t')
gene_list.to_csv(OUTPUT_GENE_LIST, index=False, sep='\t')
