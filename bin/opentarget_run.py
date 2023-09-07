#!/usr/bin/env python
# opentarget_run.py
"""""""""""""""""""""""""""""""""""""""""""""""""""""
# Author: Sarah Fong and Keila Velazquez-Arcelay
#
# unit test - Open Target Genetics PheWAS
# 
# Query open target genetics pheWAS for variant-associated traits. 
# using opentarget.py script. 
# 
# BEFORE YOU START: 
#     CHANGE OPEN_TARGET_PYPATH to str w/ full path to opentarget.py module. 
#
"""""""""""""""""""""""""""""""""""""""""""""""""""""


INPUT_FILE = '../data/opentargets_input_browning18_hg38.txt.xz'
OUTPUT_FILE = '../data/opentargets_output_introgressed_variants.tsv.xz'


import pandas as pd
import os, sys
import lzma
import opentarget as otg 
OPEN_TARGET_PYPATH="./" #<<full_path_to_opentarget.py>>
sys.path.append(os.path.join(OPEN_TARGET_PYPATH, 'opentarget.py')) 


# List of lookup variants. Format: CHR_POS_REF_ALT, str
snp_list = lzma.open(INPUT_FILE, 'rt', encoding='utf-8').readlines()[1:]
snp_list = [n.strip() for n in snp_list]
#snp_list = ["1_154453788_C_T", "16_23192369_C_T"]

# Running the OTG PheWAS function returns a Pandas dataframe containing associations. 
results = otg.run_otg_phewas(snp_list) 


# SAVE DATA
if os.path.exists(OUTPUT_FILE):
    print('FILE NOT SAVED.\n    FILENAME ALREADY EXISTS IN THIS DIRECTORY.')
    pass
else:
    results.to_csv(OUTPUT_FILE, sep='\t', index=False, compression='xz')
#results.to_csv(OUTPUT_FILE, sep='\t', index=False)

