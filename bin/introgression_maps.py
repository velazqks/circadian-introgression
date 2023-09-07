#!/usr/bin/env python
# introgression_maps.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# Date: 2023/05
#
# Description: Creates a dataset of all the available archaic introgression maps.
#              Study that also merged introgression maps: 
#              https://www.science.org/doi/10.1126/sciadv.abc0776
# 
# Introgression maps:
#  - Browning & Boeckx, 2018
#  - Sankararaman et al., 2014
#  - Schaefer et al., 2021. SARGE
#  - Skov et al., 2020. hg38. Lifted to hg19. 
#  - Steinrücken et al., 2018
#  - Vernot et al., 2016
#  - Hubisz et al., 2020. They validate their method only on a few individuals
#
"""""""""""""""""""""""""""""""""


# INPUT DATA FILES
data_files = {
    'browning2018': '../data/raw_introgressed_browning2018.bed.xz',
    'sankararaman2014': '../data/raw_introgressed_sankararaman2014.bed.xz',
    'schaefer2021': '../data/raw_introgressed_schaefer2021.bed.xz',
    'skov2020': '../data/raw_introgressed_skov2020.bed.xz',
    'steinruecken2018': '../data/raw_introgressed_steinruecken2018.bed.xz',
    'vernot2016': '../data/raw_introgressed_vernot2016.bed.xz'
}

# OUTPUT DATA
OUTPUT_FILE = '../data/introgression_maps.bed.xz'


import os,io,re
import pandas as pd
import pybedtools
#from functools import reduce


def collapse_columns(df,col):
    cols = df.columns.values.tolist()
    cols.remove(col)
    df = df.groupby(cols)[col].apply('|'.join).reset_index()
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


# LOAD FILES INTO A DICTIONARY
dfs_d = {}
for name, filename in data_files.items():
    dfs_d[name] = pd.read_csv(filename, sep='\t', compression='xz')


# print(dfs_d.keys())


# Filter by snptype, keep first 3 columns
dfs_d['skov2020'] = dfs_d['skov2020'][dfs_d['skov2020']['snptype'].isin(['DAVsnp','linkedDAVsnp'])].iloc[:,:3].drop_duplicates()

dfs_d['browning2018'] = dfs_d['browning2018'].iloc[:,:5].drop_duplicates()
dfs_d['schaefer2021'] = dfs_d['schaefer2021'].iloc[:,:3].drop_duplicates()
dfs_d['steinruecken2018'] = dfs_d['steinruecken2018'].iloc[:,:3].drop_duplicates()


# Map of name and dataframes
dfNames_d = {'Browning2018': dfs_d['browning2018'], 
             'Sankararaman2014': dfs_d['sankararaman2014'], 
             'Schaefer2021': dfs_d['schaefer2021'], 
             'Skov2020': dfs_d['skov2020'], 
             'Steinruecken2018': dfs_d['steinruecken2018'], 
             'Vernot2016': dfs_d['vernot2016']
            }


# print(dfNames_d.keys())


intersection_d = {}
for i in dfNames_d.keys():
    intersection_d[i] = intersect_a_and_b(dfs_d['browning2018'],dfNames_d[i])
    intersection_d[i].columns = intersection_d[i].columns.str.replace('_wa', '')
    intersection_d[i] = intersection_d[i].iloc[:,:4].drop_duplicates()
    intersection_d[i][i] = 1


# print(intersection_d.keys())


# Create map of all the introgression maps
map_df = None
for df in intersection_d.values():
    if map_df is None:    # Add the first df to the empty map_df
        map_df = df
    else:
        map_df = pd.merge(map_df, df, on=['Chr','Start','End','ID'], how='outer')

# Place 0s in place of NaNs
map_df.iloc[:,4:] = map_df.iloc[:,4:].fillna(0).astype(int)

# Filter rows by the number of intersecting methods
#map_df[map_df.iloc[:,4:].isnull().sum(axis=1)==0]
#map_df[map_df.iloc[:,4:].sum(axis=1)==6]

# How many variant sites are unique to Browning2018?
#len(map_df[(map_df.iloc[:,4:5].sum(axis=1)==1) & (map_df.iloc[:,5:].sum(axis=1)==0)])


# SAVE DATA
map_df.to_csv(OUTPUT_FILE, sep='\t', index=False, compression='xz')
