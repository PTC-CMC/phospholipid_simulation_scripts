import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import pdb
import json

df = pd.read_csv('roughness_summary.csv')
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['axes.labelsize'] = 24
matplotlib.rcParams['legend.fontsize'] = 18


###########
### Roughnes ###
###########
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing['msr_mean'].values)
    ffayerr_vals.append(thing['msr_std'].values)
    #unit = unit_to_str[thing['APL_unit'].values[0]]

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing['msr_mean'].values)
    ohyerr_vals.append(thing['msr_std'].values)

bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
#ax.set_ylim([28,36])
ax.set_ylabel("MSR ({})".format('nm'))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("msr.png")

