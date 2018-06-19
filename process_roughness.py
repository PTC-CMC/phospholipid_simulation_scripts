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


df = pd.DataFrame(all_dicts, columns=['DSPC', 'oh12', 'oh16', 'oh24', 'ffa12', 
                                        'ffa16', 'ffa24', 'msr_mean', 'msr_std'])
df.to_csv('roughness_summary.csv')
