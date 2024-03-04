#!/usr/bin/env python

import subprocess
import pandas as pd
import numpy as np
import math
import sys

result = sys.argv[1]
factor = sys.argv[2]
con = sys.argv[3]
dis = sys.argv[4]
ht = int(con) + int(dis)
nht = sys.argv[5]
fct = sys.argv[6]
p_value = float(sys.argv[7])

# Preprocessing a table from PERFECTOS-APE
with open(result, 'r') as file:
    lines = file.readlines()

lines[0] = lines[0].replace('# ', '')

with open(result, 'w') as file:
    file.writelines(lines)

# Combining the table from PERFECTOS-APE and the TF table from ADASTRA
result = pd.read_csv(result, sep = '\t')
result = result[['SNP name', 'motif', 'P-value 1', 'P-value 2', 'Fold change']]

tf = pd.read_csv(factor, sep = '\t')
tf = tf[tf['ID'].isin(result['SNP name'])]
tf = tf[['ID', 'fdrp_bh_ref', 'fdrp_bh_alt', 'motif_fc']]

merged = pd.merge(result, tf, left_on = 'SNP name', right_on = 'ID')
merged = merged[['P-value 1', 'P-value 2', 'fdrp_bh_ref', 'fdrp_bh_alt', 'motif']]

# Counting alleles
concordance = merged[
    (merged[['P-value 1', 'P-value 2']].min(axis = 1) < p_value) &
    ((merged['fdrp_bh_ref'] - merged['fdrp_bh_alt']) * (merged['P-value 1'] - merged['P-value 2']) > 0)
]

discordance = merged[
    (merged[['P-value 1', 'P-value 2']].min(axis = 1) < p_value) &
    ((merged['fdrp_bh_ref'] - merged['fdrp_bh_alt']) * (merged['P-value 1'] - merged['P-value 2']) < 0)
]

no_hit = merged[
    merged[['P-value 1', 'P-value 2']].min(axis = 1) >= p_value
]

hit = merged[
    merged[['P-value 1', 'P-value 2']].min(axis = 1) < p_value
]

merged['number'] = merged['motif'].apply(lambda x: int(x.split('.')[2]))
c_0 = merged['number'].value_counts().get(0, 0)
c_1 = merged['number'].value_counts().get(1, 0)
c_2 = merged['number'].value_counts().get(2, 0)
c_3 = merged['number'].value_counts().get(3, 0)

print(f'{fct}\t{con}\t{concordance.shape[0]}\t{dis}\t{discordance.shape[0]}\t{ht}\t{hit.shape[0]}\t{nht}\t{no_hit.shape[0]}\t{c_0}\t{c_1}\t{c_2}\t{c_3}')
