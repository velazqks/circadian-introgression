#!/usr/bin/env python
# opentargets_generate_input.py
"""""""""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Generate opentargets.py input file to scrape associations from a list of variants.
#              Format: CHR_POS_REF_ALT, hg38.
# 
"""""""""""""""""""""""""""""""""""""""


# INPUT DATA
INPUT_FILE = '../data/raw_introgressed_browning2018.bed.xz'
START_HG38_COL = 'Start_hg38'
END_HG38_COL = 'End_hg38'
REF_ALT_COL = 'REF/ALT'

# OUTPUT DATA
OUTPUT_FILE = '../data/opentargets_input_browning18_hg38.txt.xz'


# https://pypi.org/project/liftover/
from liftover import ChainFile
import pandas as pd
import numpy as np


# LOAD FILES
input_file = pd.read_csv(INPUT_FILE, sep='\t', compression='xz')

# GET UNIQUE LIST OF LOCI
df_uniq = input_file.iloc[:,:3].drop_duplicates().reset_index(drop=True)

# LiftOver
converter = ChainFile('../data/hg19ToHg38.over.chain.gz', 'hg19', 'hg38')

l = []
for i in range(len(df_uniq)):
    chrom = df_uniq.iat[i, 0]
    pos = df_uniq.iat[i, 2]
    hg38 = converter[chrom][pos]
    if len(hg38)==0:
        hg38 = ('',np.nan,'')
    else:
        hg38 = hg38[0]
    l.append('{}, {}, {}'.format(chrom,pos,hg38[1]))

l = [n.split(', ') for n in l]

# CREATE DATAFRAME
df_map = pd.DataFrame(l,columns=['Chr','End',END_HG38_COL])
df_map = df_map[df_map[END_HG38_COL]!='nan']
df_map = df_map.astype({'Chr':str, 'End': int, END_HG38_COL: int})

df_hg38 = pd.merge(input_file,df_map,on=['Chr','End'])

df_hg38[START_HG38_COL] = df_hg38[END_HG38_COL] - 1

# FORMAT FILE
# Remove chr string in Chr column
df_hg38['Chr'] = df_hg38['Chr'].str.replace('chr','')

df_hg38[REF_ALT_COL] = df_hg38[REF_ALT_COL].str.replace('/','_')

# ADD NEW COLUMN
df_hg38['OpenTargets'] = df_hg38['Chr'] + '_' + df_hg38[END_HG38_COL].astype(str) + '_' + df_hg38[REF_ALT_COL]

opentargets = df_hg38['OpenTargets'].drop_duplicates()


# SAVE DATA
opentargets.to_csv(OUTPUT_FILE, sep='\t', index=False, compression='xz') # quoting=0, 
"""
with open(OUTPUT_FILE, 'w') as f:
    f.write('\n'.join(opentargets))
"""