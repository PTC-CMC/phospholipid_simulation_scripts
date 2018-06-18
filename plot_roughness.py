import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import pdb
import json

unit_to_str = { 'angstrom**2': '$\AA^2$',
                'angstrom' : '$\AA$',
                'nanometer' : 'nm',
                'nanometer**2': 'nm$^2$' ,
                'degree': 'degree'}
roughdf = pd.read_csv('roughness.csv')
index = json.load(open('index.txt', 'r'))

# Make a nice pandas dataframe that aggregates sim data
all_dicts = []
all_compositions = [{'DSPC':32, 'ffa12':32}, {'DSPC':32, 'ffa16':32}, 
        {'DSPC':32, 'ffa24':32}, {'DSPC':32, 'oh12':32}, {'DSPC':32, 'oh16':32},
        {'DSPC':32, 'oh24':32}]
for composition in all_compositions:
    composition_data = {}
    for component, number in composition.items():
        composition_data[component] = number
    composition_data['msr_mean'] = []
    composition_data['msr_std'] = []
    for i, name in enumerate(index.keys()):
        if index[name]['components'] == composition:
            sub_df = roughdf.loc[roughdf['name'] == name]
            composition_data['msr_mean'].extend(sub_df['MSR_mean'].values)
            composition_data['msr_std'].extend(sub_df['MSR_std'].values)
    from uncertainties import ufloat
    from uncertainties.umath import *

    ufloats = []
    for val, std in zip(composition_data['msr_mean'], composition_data['msr_std']):
        ufloats.append(ufloat(val, std))
    foo = np.mean(ufloats)

    #composition_data['msr_mean'] = np.mean(composition_data['msr_mean'])
    #composition_data['msr_std'] = (1/3)*np.sqrt(sum(err**2 for err in composition_data['msr_std']))
    #print("By hand: {} {}".format(composition_data['msr_mean'], composition_data['msr_std']))
    #print("Uncertainties: {}".format(foo))
    composition_data['msr_mean'] = foo.n
    composition_data['msr_std'] = foo.s

    all_dicts.append(composition_data)


df = pd.DataFrame(all_dicts)
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
    #unit = unit_to_str[thing['APL_unit'].values[0]]

bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
#ax.set_ylim([28,36])
ax.set_ylabel("MSR ({})".format('nm$^2$'))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("msr.png")

