import os
import pdb
import json
import numpy as np
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd

#####################
### From the dipole relaxation fits, average and error
##############

index = json.load(open('index.txt', 'r'))
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()

all_compositions = [{'DSPC':32, 'ffa12':32}, {'DSPC':32, 'ffa16':32}, 
        {'DSPC':32, 'ffa24':32}, {'DSPC':32, 'oh12':32}, {'DSPC':32, 'oh16':32},
        {'DSPC':32, 'oh24':32}]

for composition_i in all_compositions:
    # Initialize data structures
    composition_data = {}
    composition_data['dip_tau_mean'] = []

    composition_data['dip_tau_std'] = []

    composition_data['unit'] = 'ps'

    for component, number in composition_i.items():
        composition_data[component] = number
        for i, name in enumerate(index.keys()):
            if index[name]['components'] == composition_i:
                os.chdir(os.path.join(curr_dir, name))
                #dip_corr = np.loadtxt('dip_corr.dat')
                dip_corr = json.load(open('dip_corr_fit.txt', 'r'))
                composition_data['dip_tau_mean'].append(dip_corr['tau'])

    
    composition_data['dip_tau_std'] = np.std(composition_data['dip_tau_mean']) 
    composition_data['dip_tau_mean'] = np.mean(composition_data['dip_tau_mean']) 

    summary_df = summary_df.append(composition_data, ignore_index=True)


os.chdir(curr_dir)
summary_df.to_csv('dipole_relaxation.csv')
