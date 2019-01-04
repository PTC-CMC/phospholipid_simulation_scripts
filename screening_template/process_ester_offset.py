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

index = json.load(open('combined_index.json','r'))
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()
# Make the composition dictionary

all_compositions = [{'DSPC':21, 'oh18':43}, {'DSPC':21, 'oh20':43}, 
        {'DSPC':21, 'oh22':43}]


for composition_i in all_compositions:
    # Initialize data structures
    composition_data = OrderedDict()
    composition_data['ester_offset_mean'] = []
    composition_data['ester_offset_std'] = []
    composition_data['ester_offset_unit'] = []

    # Aggregate all the associated raw data from each simulation
    for i, name in enumerate(index.keys()):
        os.chdir(os.path.join(curr_dir, name))
        sim_data = json.load(open('ester_offset.txt','r'), object_pairs_hook=OrderedDict)

        if index[name]['components'] == composition_i:
            composition_data['ester_offset_mean'].append(sim_data['ester_offset_mean'])
            composition_data['ester_offset_std'].append(sim_data['ester_offset_std'])
            composition_data['ester_offset_unit'] = sim_data['ester_offset_unit']
            for component, number in composition_i.items():
                composition_data[component] = number


    
    all_df = all_df.append(composition_data, ignore_index=True)

    # Now compute averages and propagate error
    ufloats = [ufloat(val,std) for val,std in 
            zip(composition_data['ester_offset_mean'], 
                  composition_data['ester_offset_std'])]
    foo = np.mean(ufloats)
    composition_data['ester_offset_mean'] = foo.n
    composition_data['ester_offset_std'] = foo.s


    summary_df = summary_df.append(composition_data, ignore_index=True)

os.chdir(curr_dir)
summary_df.to_csv('ester_offset.csv')

