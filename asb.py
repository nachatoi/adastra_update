#!/usr/bin/env python

import pandas as pd
import sys

factor = str(sys.argv[1])
p = float(sys.argv[2])
short_name = str(sys.argv[3])

tf = pd.read_csv(factor, sep = '\t')
tf = tf[tf[['fdrp_bh_ref', 'fdrp_bh_alt']].min(axis = 1) < p]
tf.to_csv('filtered/' + short_name + '.filtered', sep = '\t', index = False)
