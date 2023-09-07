#!/usr/bin/env python
# input_for_morningness_cumulative_fraction.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Create cumulative fraction of p-values associated with the UK Biobank 
#              morningness phenotyme in a GWAS analysis.
#              
# 
"""""""""""""""""""""""""""


# INPUT DATA
MORNINGNESS_FILE = '../data/raw_nealelab_round2_1180.bed.xz'
BROWNING_FILE = '../data/introgression_maps.bed.xz'
CIRCADIAN_INTROGRESSED = '../data/circadian_variants_introgressed.bed.xz'


import pandas as pd
import numpy as np
import math
import re
import subprocess
import tempfile
import warnings
warnings.filterwarnings('ignore', category=pd.core.common.SettingWithCopyWarning)



def run_all_clump(threshold,ld,out,out_c):
    # Extract introgressed chronotype SNPs
    data = pd.merge(morningness, maps, on=['Chr','Start','End'])
    
    # Filter introgressed SNPs by level of support of detection methods
    data = data[data.iloc[:,-6:].sum(axis=1)>=threshold].iloc[:,:-6]
    
    # Get a subset of SNPs from circadian genes
    data_c = pd.merge(introgressed_circadian, data, on=['Chr','Start','End'])
    
    # LD clump both sets
    dfs = {out:data, out_c:data_c}
    df_l = []
    for k,v in dfs.items():
        # Create the Plink Clump input
        to_clump = create_clump_input(v)
        
        # Clump sets of variants by R2>=threshold
        clump(to_clump, ld, k)

        # Wrangle clump results
        clumped = wrangle_clumped(v, k)

        # Cumulative fraction of morningness effect
        cumulative_frac = cumulative_fraction(clumped)
        
        # Append each df to a list
        df_l.append(cumulative_frac)
        
    return df_l

def create_clump_input(data):
    df = data[['ID','pval']]
    df.rename(columns={'ID':'SNP','pval':'P'}, inplace=True)
    return df

def clump(df,ld,output):
    # Write the data to a temporary file
    temp_dir = '../data/'
    temp_file = tempfile.NamedTemporaryFile(suffix='.txt', delete=True, dir=temp_dir)
    print(temp_file.name)
    
    df.to_csv(temp_file.name, sep='\t', index=False)

    # Plink will define the lead SNPs by the highest P-value in a haplotype of LD=0.5, 
    # and the SNPs in LD with the lead SNP are secondary or proxy SNPs
    plink_cmd = [
        'plink',
        '--bfile',
        '../general_data/g1000_eur_introgressed',
        '--clump',
        temp_file.name,
        '--clump-r2',
        '0.50',
        '--clump-field',
        'P',
        '--clump-p1',
        '1',
        '--out',
        temp_dir+output
    ]
    
    # Execute the Plink command and capture the output
    # NOTE: My HPC is running on Python 3.6.11, and 'capture_output' isn't available
    # result = subprocess.run(plink_cmd, check=True, capture_output=True, text=True)
    result = subprocess.run(plink_cmd, check=True, stdout=subprocess.PIPE)
    
    # Get the output as a string
    output_bytes = result.stdout
    output_string = output_bytes.decode()

    # Remove the temporary file
    temp_file.close()

    
def wrangle_clumped(df, file):
    # Import clumped output
    lst = [n.strip() for n in open(f'../data/{file}.clumped').readlines()]
    lst = [re.split(r'\s+', m) for m in lst]#{2,}
    # Remove empty lists at the end
    lst = [l for l in lst if not all(element == '' for element in l)]
    
    # Create a dataframe from the list of chronotype SNPs
    data = pd.DataFrame(lst[1:], columns=lst[0])
    data.rename(columns={'SNP':'ID'}, inplace=True)
    
    # Sort the df rows by chromosome and position
    data[['CHR','BP']] = data[['CHR','BP']].astype(int)
    data.sort_values(by=['CHR','BP'], inplace=True)
    data.reset_index(drop=True, inplace=True)
    data = data[['ID','P','SP2']].drop_duplicates()
    data['P'] = data['P'].astype(float)
    
    # Merge to add beta values
    data = pd.merge(df, data, on=['ID'])
    
    return data


def cumulative_fraction(data):
    # Cumulative fraction of effects on morningness by SNP significance
    df = data.sort_values('P')
    df=df.reset_index(drop=True)
    c=0    # Number of chronotype SNPs with positive beta value
    l=[]
    for i in range(len(df)):
        if df['beta'][i] >= 0:
            c+=1
        l.append([df['ID'][i],
                  df['beta'][i],
                  (df['P'][i]),
                  math.log10(df['P'].astype(float)[i]),
                  c/(i+1)])
    df_to_plot = pd.DataFrame(l,columns=['RSID','Beta','P-Value','P-Value(log10)','Cumulative_Fraction'])
    return df_to_plot


# LOAD FILES
# Load UKBiobank Morningness SNPs
morningness = pd.read_csv(MORNINGNESS_FILE, sep='\t').iloc[:,:-1]

# Check for duplicated loci in the morningness df
#morningness[morningness[['Chr','Start','End']].duplicated(keep=False)]

# Drop duplicates based on Locus and keep the lowest pval
morningness = morningness.sort_values(['Chr','Start','End','pval'])
morningness = morningness.drop_duplicates(subset=['Chr', 'Start', 'End'], keep='first')
morningness = morningness.reset_index(drop=True)

# Load introgression maps table
maps = pd.read_csv(BROWNING_FILE, sep='\t')
maps = maps[maps['ID'].str.startswith('rs')]

# Introgressed SNPs in circadian genes
introgressed_circadian = pd.read_csv(CIRCADIAN_INTROGRESSED, sep='\t').iloc[:,:3].drop_duplicates()


# INTROGRESSED VARIANTS IDENTIFIED BY BROWNING ET AL., 2018 AND AT LEAST ONE OTHER METHOD
to_plot_browning = run_all_clump(1, 0.5, 
                                  'introgressed_morningness_browning2018',
                                  'introgressed_morningness_c_browning2018')

# INTROGRESSED VARIANTS IDENTIFIED BY BROWNING ET AL., 2018 AND AT LEAST ONE OTHER METHOD
to_plot_2_methods = run_all_clump(2, 0.5, 
                                  'introgressed_morningness_support_by_2',
                                  'introgressed_morningness_c_support_by_2')

# INTROGRESSED VARIANTS IDENTIFIED BY ALL THE INTROGRESSION DETECTION METHODS
to_plot_all_methods = run_all_clump(6, 0.5, 
                                  'introgressed_morningness_support_by_all',
                                  'introgressed_morningness_c_support_by_all')


# SAVE DATA
dataframes = {
    '../data/plotting_morningness_cumulative_fraction_browning.tsv':to_plot_browning[0],
    '../data/plotting_morningness_cumulative_fraction_browningc.tsv':to_plot_browning[1],
    '../data/plotting_morningness_cumulative_fraction_i2.tsv':to_plot_2_methods[0],
    '../data/plotting_morningness_cumulative_fraction_i2c.tsv':to_plot_2_methods[1],
    '../data/plotting_morningness_cumulative_fraction_iall.tsv':to_plot_all_methods[0],
    '../data/plotting_morningness_cumulative_fraction_iallc.tsv':to_plot_all_methods[1]
}

for filename,df in dataframes.items():
    df.to_csv(filename, sep='\t', index=False)
