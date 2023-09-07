#!/usr/bin/env python
# table_circadian_genes.py
"""""""""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Webscrape descriptions for each of the candidate circadian genes.
#              Create Table 1 containing circadian genes, descriptions, and evidence sources.
#
"""""""""""""""""""""""""""""""""


# INPUT DATA
CIRCADIAN_GENES = '../data/circadian_genes.list'
CANDIDATE_GENES = '../data/circadian_genes_candidate.tsv'


# OUTPUT DATA
OUTPUT_FILE = '../data/table_circadian_genes.tsv'


import pandas as pd
import requests, sys


def transform_evidence_source_columns(row):
    #return [col for col in row.index if row[col] == 1]
    return ', '.join(row.index[row == 1])

def scrape_ensembl(genes,data):
    # Scrape gene description
    # https://rest.ensembl.org/documentation/info/xref_name

    l = []
    for gene in genes:
        server = 'https://rest.ensembl.org'
        ext = '/xrefs/id/{}?'.format(gene)

        r = requests.get(server+ext, headers={'Content-Type':'application/json'})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        l.append('{}; {}'.format(
            (gene),
            (decoded[3]['description'])
        ))
    
    # Split list values
    l = [n.split('; ') for n in l]
    # Create df
    df = pd.DataFrame(l, columns=['GeneID', 'Description'])
    # Merge
    df = pd.merge(data, df, on='GeneID').sort_values(by='GeneID')
    
    return df


# IMPORT CIRCADIAN GENES LIST
circadian_genes = pd.read_csv(CIRCADIAN_GENES, sep='\t')
candidate_genes = pd.read_csv(CANDIDATE_GENES, sep='\t')
candidate_genes = candidate_genes[candidate_genes['Confidence']!='Low']

# CREATE COLUMN SUMMARIZING THE SOURCES OF EVIDENCE FOR EACH GENE
candidate_genes['Evidence'] = candidate_genes.iloc[:, 3:].apply(transform_evidence_source_columns, axis=1)
#candidate_genes['Evidence'] = candidate_genes.iloc[:,3:].stack().groupby(level=0).agg(', '.join)
candidate_genes = candidate_genes[['GeneID', 'GeneName', 'Evidence', 'Confidence']]

# WEBSCRAPE GENE DESCRIPTIONS FOR EACH OF THE CIRCADIAN GENES
table = scrape_ensembl(circadian_genes['GeneID'].values.tolist(),
                                 circadian_genes)

# MERGE GENE DESCRIPTIONS WITH EVIDENCE SOURCES AND CONFIDENCE LEVEL
table = pd.merge(table, candidate_genes, on=['GeneID','GeneName']).sort_values(by='GeneName')


# SAVE DATA
table.to_csv(OUTPUT_FILE, sep='\t', index=False)
#table.to_excel('../data/table_1.xlsx', index=False)
