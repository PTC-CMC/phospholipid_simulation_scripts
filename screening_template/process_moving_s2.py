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
### process trials of moving s2 windows
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
        # this dictionary should keep track of the nubmer of each component
        composition_data[component] = number
        # Look through each simulation to try and find a match
    for i, name in enumerate(index.keys()):
        if index[name]['components'] == composition_i:
            os.chdir(os.path.join(curr_dir, name))
            moving_s2 = json.load(open('moving_s2.json','r'))
            for key in composition_i.keys():
                if 'DSPC' == key:
                    if (key + "a_s2") not in composition_data.keys():
                        composition_data[key+"a_s2"] = []
                    if (key + "b_s2") not in composition_data.keys():
                        composition_data[key+"b_s2"] = []

                    composition_data[key+"a_s2"].append(
                            moving_s2[key+"a"]['s2_mean'])
                    composition_data[key+"b_s2"].append(
                            moving_s2[key+"b"]['s2_mean'])
                else:
                    if (key+'_s2') not in composition_data.keys():
                        composition_data[key+'_s2'] = []

                    composition_data[key+'_s2'].append(
                        moving_s2[key]['s2_mean'])

    # Compute averages and errors
    s2_keys = [key for key in composition_data.keys() if '_s2' in key]
    for key in s2_keys:
        arr = np.array(composition_data[key])
        mean = np.mean(arr,axis=0)
        err = (np.std(arr,axis=0) 
                / np.sqrt(np.shape(arr)[0]))
        err_key = key[:-2]+"err"
        composition_data[key] = mean
        composition_data[err_key] = err
    summary_df = summary_df.append(composition_data, ignore_index=True)
os.chdir(curr_dir)
summary_df.to_csv('moving_s2.csv')
