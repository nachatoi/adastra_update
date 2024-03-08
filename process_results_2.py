#!/usr/bin/env python

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import warnings
import shutil
import sys 
import os

p_value = float(sys.argv[1])
adastra = sys.argv[2]

warnings.simplefilter(action = 'ignore', category = Warning)

def concordance(result, factor, con, dis, nht, fct, p_value):
    p_value = float(p_value)
    ht = int(con) + int(dis)

    # Preprocessing a table from PERFECTOS-APE
    with open(result, 'r') as file:
        lines = file.readlines()

    lines[0] = lines[0].replace('# ', '')

    with open(result, 'w') as file:
        file.writelines(lines)

    # Combining the table from PERFECTOS-APE and the TF table from ADASTRA
    result = pd.read_csv(result, sep = '\t')
    result = result[['SNP name', 'P-value 1', 'P-value 2', 'min P-value', 'Fold change', 'motif', 'Concordance', 'Discordance']]

    tf = pd.read_csv(factor, sep = '\t')
    tf = tf[tf['ID'].isin(result['SNP name'])]
    tf = tf[['ID', 'fdrp_bh_ref', 'fdrp_bh_alt']]

    merged = pd.merge(result, tf, left_on = 'SNP name', right_on = 'ID')

    for index, row in merged.iterrows():
        if merged.at[index, 'min P-value'] < p_value:
            smth = ((merged.at[index, 'fdrp_bh_ref'] - merged.at[index, 'fdrp_bh_alt']) * (merged.at[index, 'P-value 1'] - merged.at[index, 'P-value 2']))
            if smth > 0: 
                merged.at[index, 'Concordance'] = True
            if smth < 0: 
                merged.at[index, 'Discordance'] = True

    merged = merged.sort_values(by = ['Concordance'])
    merged.drop_duplicates(subset = 'SNP name', keep = 'last', inplace = True)

    # Counting alleles
    concordance = merged[merged['Concordance'] != False]
    discordance = merged[merged['Discordance'] != False]
    no_hit = merged[(merged['Concordance'] == False) & (merged['Discordance'] == False)]
    hit = merged[(merged['Concordance'] != False) | (merged['Discordance'] != False)]

    c_0 = merged['motif'].value_counts().get(0, 0)
    c_1 = merged['motif'].value_counts().get(1, 0)
    c_2 = merged['motif'].value_counts().get(2, 0)
    c_3 = merged['motif'].value_counts().get(3, 0)

    return f'{fct}\t{con}\t{concordance.shape[0]}\t{dis}\t{discordance.shape[0]}\t{ht}\t{hit.shape[0]}\t{nht}\t{no_hit.shape[0]}\t{c_0}\t{c_1}\t{c_2}\t{c_3}'

folders = ['pwm_results_0', 'pwm_results_1', 'pwm_results_2', 'pwm_results_3']
file_paths = dict()

for folder in folders:
    files = os.listdir(folder)
    for file in files:
        file_path = os.path.join(folder, file)
        if file in file_paths: 
            file_paths[file].append(file_path)
        else: 
            file_paths[file] = [file_path]

output_folder = 'pwm_results'
os.makedirs(output_folder, exist_ok = True)

for file, paths in file_paths.items():
    df = pd.DataFrame()
    for path in paths: 
        index = path.split('/')
        index = index[0].split('_')
        index = index[-1]
        data = pd.read_csv(path, sep = '\t')
        data['motif'] = index
        data['Concordance'] = False
        data['Discordance'] = False
        df = pd.concat([df, data], ignore_index = True)
    df['min P-value'] = df[['P-value 1', 'P-value 2']].min(axis = 1)

    file_name = os.path.basename(file)
    dest_path = os.path.join(output_folder, file_name)
    df.to_csv(dest_path, sep = '\t', index = False)

with open('adastra_pwm.tsv', mode = 'w') as pwm: 
    pwm.write('Factor\tAdastra_concordance\tConcordance\tAdastra_discordance\tDiscordance\tAdastra_hit\tHit\tAdastra_ho_hit	No_hit\t0\t1\t2\t3\n')

files = os.listdir(output_folder)

for file in files: 
    factor = file.split('.')
    factor = factor[0]
    file = output_folder + '/' + file

    for filename in os.listdir(adastra):
        if filename.startswith(f'{factor}_'):
            name = adastra + filename

            data = pd.read_csv(f'filtered/{factor}.filtered', sep = '\t')
            data['motif_conc'] = data['motif_conc'].astype(str).fillna('')

            conc = int(data['motif_conc'].str.count('Concordant').sum())
            disc = int(data['motif_conc'].str.count('Discordant').sum())
            nohit = int(data['motif_conc'].str.count('No Hit').sum())

            output = concordance(file, name, str(conc), str(disc), str(nohit), factor, str(p_value))

            with open('adastra_pwm.tsv', mode = 'a') as pwm: 
                pwm.write(output + '\n')
