import glob
import os
from collections import OrderedDict
import pdb
import json
import numpy as np
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd

#####################
### process trials of lipid rotations 
##############

index = json.load(open('new_index.txt', 'r'))
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()

all_compositions = []
for oh in ['oh12', 'oh16', 'oh24']:
    for ffa in ['ffa12', 'ffa16', 'ffa24']:
        comp = OrderedDict()
        comp['DSPC'] = 22
        comp[oh] = 21
        comp[ffa] = 21
        all_compositions.append(comp)



for composition_i in all_compositions:
    composition_data = {}
    # Initialize data structures
    for component, number in composition_i.items():
        composition_data[component] = number
        for i, name in enumerate(index.keys()):
            if index[name]['components'] == composition_i:
                os.chdir(os.path.join(curr_dir, name))
                all_scd_files = glob.glob('*_scd.xvg')
                # For each scd.xvg file, need to initialize the dictionary entry
                for scd in all_scd_files:
                    data = np.loadtxt(scd, comments=['@', '#'])
                    key = scd[:-4] + "_mean"
                    if key not in composition_data.keys():
                        composition_data[key] = []
                    composition_data[key].append(data[:,1])
    # Compute averages and errors
    scd_keys = [key for key in composition_data.keys() if 'scd_mean' in key]
    for key in scd_keys:
        if 'scd' in key:
            arr = np.array(composition_data[key])
            mean = np.mean(arr,axis=0)
            err = (np.std(arr,axis=0) 
                    / np.sqrt(np.shape(arr)[0]))
            err_key = key[:-4]+"err"
            composition_data[key] = mean
            composition_data[err_key] = err
    summary_df = summary_df.append(composition_data, ignore_index=True)
os.chdir(curr_dir)
summary_df.to_csv('scd_order.csv')
