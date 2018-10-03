import os
import json
import numpy as np
from collections import OrderedDict
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd
import pdb
################
## Script to aggregate simulation data by composition
################

index = json.load(open('index.txt','r'))
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()
# Make the composition dictionary
all_compositions = []
for oh in ['oh12', 'oh16', 'oh24']:
    for ffa in ['ffa12', 'ffa16', 'ffa24']:
        comp = OrderedDict()
        comp['DSPC'] = 22
        comp[oh] = 21
        comp[ffa] = 21
        all_compositions.append(comp)

#all_compositions = [{'DSPC':32, 'ffa12':32}, {'DSPC':32, 'ffa16':32}, 
#        {'DSPC':32, 'ffa24':32}, {'DSPC':32, 'oh12':32}, {'DSPC':32, 'oh16':32},
#        {'DSPC':32, 'oh24':32}]
for composition_i in all_compositions:
    # Initialize data structures
    composition_data = OrderedDict()
    composition_data['interface_idig_mean'] = []
    composition_data['interface_idig_std'] = []
    composition_data['interface_idig_unit'] = []

    # Aggregate all the associated raw data from each simulation
    for i, name in enumerate(index.keys()):
        os.chdir(os.path.join(curr_dir, name))
        sim_data = json.load(open('interface_idig.txt','r'), object_pairs_hook=OrderedDict)

        if index[name]['components'] == composition_i:
            composition_data['interface_idig_mean'].append(sim_data['interface_idig_mean'])
            composition_data['interface_idig_std'].append(sim_data['interface_idig_std'])
            composition_data['interface_idig_unit'] = sim_data['interface_idig_unit']
            for component, number in composition_i.items():
                composition_data[component] = number


    
    all_df = all_df.append(composition_data, ignore_index=True)

    # Now compute averages and propagate error
    ufloats = [ufloat(val,std) for val,std in 
            zip(composition_data['interface_idig_mean'], 
                  composition_data['interface_idig_std'])]
    foo = np.mean(ufloats)
    composition_data['interface_idig_mean'] = foo.n
    composition_data['interface_idig_std'] = foo.s


    summary_df = summary_df.append(composition_data, ignore_index=True)

os.chdir(curr_dir)
summary_df.to_csv('interface_idig.csv')
#all_df.to_csv('all_sims.csv')
