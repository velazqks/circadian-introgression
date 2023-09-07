#!/usr/bin/env python
# parse_kuhlwilm2019.py
"""""""""""""""""""""""""""""""""
# Author: John A. Capra
# Date: 2022/03/16
# 
# Description: Parse the human and archaic variant sites from Kuhlwilm and Boeckx 2019. 
#              Writes tab separated files with human and archaic specific variants
#              Note that these don't exactly match the counts from the manuscript.
#              Run: parse_kuhlwilm2019.py ../data/raw_kuhlwilm2019_filtered.bed.xz
#
# Source: "A catalog of single nucleotide changes distinguishing modern humans from archaic hominins"
#         https://www.nature.com/articles/s41598-019-44877-x
#         https://figshare.com/articles/dataset/Variants_and_annotations_of_Neandertals/8184038
#
"""""""""""""""""""""""""""""""""

import os
import sys
import csv
import lzma

header = None
fixed_archaic = []
fixed_human = []
fixed_human_tony = [] # preferred stricter definition

# with open(sys.argv[1], 'rb') as csvfile:
with lzma.open(sys.argv[1], 'rt') as csvfile:
    r = csv.reader(csvfile, delimiter='\t')
    for row in r:
        row = [float('NaN') if x == '' else x for x in row]
        if row[3] == 'POS':
            header = row
            continue

        pos = row[3]
        
        if row[5] == 'NaN': continue
        human_daf = float(row[5])
        
        ref = row[6]
        alt = row[7]
        if alt not in ['A','T','C','G']: continue # consider only bi-allelic

        # consider only sites with a confident ancestral allele
        ancestral_allele = row[8]
        if ancestral_allele not in ['A','T','C','G']: continue
        derived_allele = alt if ancestral_allele == ref else ref

        altai_gt = row[9]
        altai_filter = int(row[10])
        altai_allele = row[11] if altai_filter else ''

        vindija_gt = row[12]
        vindija_filter = int(row[13])
        vindija_allele = row[14] if vindija_filter else ''

        denisova_gt = row[15]
        denisova_filter = int(row[16])
        denisova_allele = row[17] if denisova_filter else ''

        arch_filt_derived_count = float(row[18])

        # define archaic specific
        if human_daf < 0.00001:
            if arch_filt_derived_count == 3:
                fixed_archaic.append(row)

        
        # define human specific
        if human_daf == 1.0:
                # require at least two archaics to pass the filter
                num_archaic_filt = altai_filter + vindija_filter + denisova_filter
                if num_archaic_filt < 2: continue

                # require all archaics to be ancestral
                arch_anc_count = 0.
                if altai_filter and altai_allele == ancestral_allele:
                    arch_anc_count += 1
                if vindija_filter and vindija_allele == ancestral_allele:
                    arch_anc_count += 1
                if denisova_filter and denisova_allele == ancestral_allele:
                    arch_anc_count += 1

                if arch_anc_count / num_archaic_filt == 1:
                    fixed_human.append(row)

        # define human specific tony-style
        if human_daf == 1.0:
                # require all three archaics to pass the filter
                num_archaic_filt = altai_filter + vindija_filter + denisova_filter
                if num_archaic_filt < 3: continue

                # require all archaics to be homozygous for ancestral
                arch_anc_count = 0.
                if altai_filter and altai_gt in ['0', '1'] and altai_allele == ancestral_allele:
                    arch_anc_count += 1
                if vindija_filter and vindija_gt in ['0', '1'] and vindija_allele == ancestral_allele:
                    arch_anc_count += 1
                if denisova_filter and denisova_gt in ['0', '1'] and denisova_allele == ancestral_allele:
                    arch_anc_count += 1

                if arch_anc_count / num_archaic_filt == 1:
                    fixed_human_tony.append(row)


print("    Num. fixed archaic:", len(fixed_archaic))
print("      Num. fixed human:", len(fixed_human))
print("Num. fixed human, Tony:", len(fixed_human_tony))


# Write the output files
with open('../data/raw_kuhlwilm19_archaic_fixed.bed', 'w') as csvfile:
    w = csv.writer(csvfile, delimiter='\t')
    w.writerow(header)
    for row in fixed_archaic: w.writerow(row)

with open('../data/raw_kuhlwilm19_human_fixed.bed', 'w') as csvfile:
    w = csv.writer(csvfile, delimiter='\t')
    w.writerow(header)
    for row in fixed_human: w.writerow(row)

with open('../data/raw_kuhlwilm19_human_fixed_tony.bed', 'w') as csvfile:
    w = csv.writer(csvfile, delimiter='\t')
    w.writerow(header)
    for row in fixed_human_tony: w.writerow(row)
