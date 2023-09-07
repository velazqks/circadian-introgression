#!/usr/bin/env python
# input_for_morningness_latitude_cline.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Extract circadian introgressed variants associated with chronotype 
#              in a UKBiobank GWAS. Population latitudes used for linear regression 
#              analysis to identify SNPs in this set following a latitudinal cline 
#              in 1000 Genomes Projects populations from Eurasian ancestry.
# 
#              # WARNING: Runs for ~45 minutes
# 
"""""""""""""""""""""""""""""""""


# INPUT DATA
KGP_META = '../data/raw_KGP_metadata.tsv'
MORNINGNESS = '../data/raw_nealelab_round2_1180.bed.xz'
INTROGRESSED = '../data/introgression_maps.bed.xz'
CHR = 'chr2'
LDLINK_TOKEN = '' # How to get LDLink token: https://ldlink.nih.gov/?tab=apiaccess

filenames = {'browning2018':'../data/introgressed_morningness_browning2018',
            'browning2018_c':'../data/introgressed_morningness_c_browning2018',
            'browning2018_plus1':'../data/introgressed_morningness_support_by_2',
            'browning2018_plus1_c':'../data/introgressed_morningness_c_support_by_2',
            'browning2018_plusAll':'../data/introgressed_morningness_support_by_all',
            'browning2018_plusAll_c':'../data/introgressed_morningness_c_support_by_all'
}


import pandas as pd
import re
import requests
import time


def wrangle_clumped(FILENAME):
    # Import ld-clumped chronotype SNPs
    with open(FILENAME+'.clumped') as file:
        l = [re.split(r'\s+', line.strip()) for line in file.readlines()]
    # Remove empty nested lists at the end
    l = [n for n in l if n != ['']]

    # Create a dataframe from the list of chronotype SNPs
    haps = pd.DataFrame(l[1:], columns=l[0])[['CHR','BP','SNP','P','SP2']]
    haps.rename(columns={'SNP':'ID'}, inplace=True)
    
    # Sort the rows by chromosome and position
    haps[['CHR','BP']] = haps[['CHR','BP']].astype(int)
    haps.sort_values(by=['CHR','BP'], inplace=True)
    haps.reset_index(drop=True, inplace=True)
    #haps = haps[['ID','SP2']].drop_duplicates()
    
    # Make the index a column to create haplotype ids from it
    haps.reset_index(drop=False, inplace=True)
    haps.rename(columns={'index': 'Haplotype'}, inplace=True)
    # Add 1 to each index number to start from 1 instead of 0
    haps['Haplotype'] = haps['Haplotype'] + 1
    
    # Find the max number of digits in the haplotype ids to match the number of digits   
    num_digits = len(str(haps['Haplotype'].max()))
    # Pad the numbers in the Haplotype column with leading zeros
    haps['Haplotype'] = haps['Haplotype'].astype(str).str.zfill(num_digits)
    
    # Explode the values in the column containing the SNPs in the haplotype
    haps['SP2'] = haps['SP2'].str.split(',')
    haps = haps.explode('SP2')
    
    # Extract the RSIDs from the exploded column and sort by chromosome and position
    haps['SP2'] = haps['SP2'].str.extract(r'(rs\d+)')
    
    # Create a df containing all the SNPs in each haplotype, haplotype id, and p-value
    lead_snps = haps[['Haplotype','ID']].drop_duplicates()
    proxy_snps = haps[['Haplotype','SP2']].drop_duplicates().dropna()
    proxy_snps.rename(columns={'SP2':'ID'}, inplace=True)
    snps = pd.concat([lead_snps, proxy_snps])
    snps.reset_index(drop=True, inplace=True)
    snps = snps.sort_values('Haplotype').reset_index(drop=True)
    
    # Merge with initial plink input file to add loci and beta
    snps = pd.merge(imorningness, snps, on='ID')
    snps.rename(columns={'ID':'RSID'}, inplace=True)
    
    return snps



def allele_frq(df, token):
    # https://ldlink.nih.gov/?tab=apiaccess
    url = 'https://ldlink.nih.gov/LDlinkRest/ldpop'
    
    not_found = []
    pop_frq_d = {}
    for i in df:
        try:
            params = {
                'var1': i,
                'var2': i,
                'pop': 'CDX+CHB+JPT+KHV+CHS+BEB+GIH+ITU+PJL+STU+GBR+FIN+IBS+TSI+CEU',
                'r2_d': 'r2',
                'genome_build': 'grch37',
                'token': token
            }
            
            response = requests.get(url, params=params)
            # print(response.content)
            
            if response.status_code == 200:
                lol = [n.split('\t') for n in response.content.decode().split('\n') if n != '']
                pop_frq = pd.DataFrame(lol[1:], columns=lol[0]).iloc[:,[1,3]]
                pop_frq['RSID'] = i
                #pop_frq[['REF_FREQ', 'ALT_FREQ']] = pop_frq.iloc[:,1].str.split(', ', expand=True)
                split_result = pop_frq.iloc[:, 1].str.split(', ', expand=True)
                pop_frq['REF_FREQ'] = split_result[0]
                pop_frq['ALT_FREQ'] = split_result[1]
                pop_frq.drop(pop_frq.columns[1], axis=1, inplace=True)
                pop_frq.rename(columns={pop_frq.columns[0]:'Population'}, inplace=True)
                pop_frq = pop_frq.replace('[A-Z]+: |%','',regex=True)
                pop_frq_d[i] = pop_frq
            else:
                print('Error:', response.status_code)
            
            print(i)
            time.sleep(1)
        except:
            print(f'{i} was not found')
            not_found.append(i)
        
    return pop_frq_d, not_found



# LOAD DATA
morningness = pd.read_csv(MORNINGNESS, sep='\t', compression='xz').iloc[:,:-1]
introgressed = pd.read_csv(INTROGRESSED, sep='\t', compression='xz')
imorningness = pd.merge(morningness, introgressed, on=['Chr','Start','End'])[
               ['Chr','Start','End','ID','beta','pval']]
lat = pd.read_csv(KGP_META, sep=',', index_col=0).iloc[:,[0,2,-1]]
lat.rename(columns={'pop':'Population'}, inplace=True)

# LD-clumped chronotype introgressed SNPs.
confidence_sets = {}
for k,v in filenames.items():
    output = wrangle_clumped(v)
    confidence_sets[k] = output


# Webscrape Allele frequencies
# WARNING: Runs for ~45 minutes
start_time = time.time()
first_item = confidence_sets[next((k) for k,v in filenames.items())]
frqs, no_frq = allele_frq(first_item[first_item['Chr']==CHR]['RSID'], LDLINK_TOKEN)
frqs_df = pd.concat(frqs.values(), ignore_index=True)
#frqs = pd.read_csv('../data/introgressed_morningness_allele_frq.tsv', sep='\t')
end_time = time.time()

print('Duration: {} minutes'.format((end_time - start_time)/60))


# Add latitude to the df 
lat_frqs = pd.merge(lat,frqs_df, on='Population')

# Update the sets with population latidudes and allele frequencies 
for k,v in confidence_sets.items():
    confidence_sets[k] = pd.merge(lat_frqs, v['RSID'], on='RSID')
    confidence_sets[k] = pd.merge(first_item,v,on='RSID')


# SAVE DATA
for k,v in confidence_sets.items(): 
    confidence_sets[k].to_csv(f'../data/plotting_morningness_latitude_cline_{k}_{CHR}.bed', 
                           sep='\t', index=False)

# Optional: Convert column data types
#result['REF_FREQ'] = pd.to_numeric(result['REF_FREQ'], errors='coerce')
#data['ALT_FREQ'] = pd.to_numeric(result['ALT_FREQ'], errors='coerce')
