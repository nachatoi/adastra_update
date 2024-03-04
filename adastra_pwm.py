#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import sys

warnings.simplefilter(action = 'ignore', category = RuntimeWarning)

output = sys.argv[1]

pwm_0 = pd.read_csv('adastra_pwm_0.tsv', sep='\t')
pwm_1 = pd.read_csv('adastra_pwm_1.tsv', sep='\t')
pwm_2 = pd.read_csv('adastra_pwm_2.tsv', sep='\t')
pwm_3 = pd.read_csv('adastra_pwm_3.tsv', sep='\t')

merged = pd.concat([pwm_0, pwm_1, pwm_2, pwm_3])

merged['Log_Adastra_Conc_Disc'] = np.log10(merged['Adastra_concordance'] / merged['Adastra_discordance'])
merged['Log_Conc_Disc'] = np.log10(merged['Concordance'] / merged['Discordance'])
merged['delta'] = merged['Log_Adastra_Conc_Disc'] - merged['Log_Conc_Disc']

merged.sort_values(by='delta', ascending=True, inplace=True)
pwm = merged.drop_duplicates(subset='Factor', keep='first')
pwm.to_csv(output, sep='\t', index=False)
