import os
import pdb
import json
import numpy as np
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd

#####################
### From the various water orientation relaxation timeseries,
### Integrate to get the respective lifetime, tau
##############

index = json.load(open('index.txt', 'r'))
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()

all_compositions = [{'DSPC':21, 'ffa12':43}, {'DSPC':21, 'ffa16':43}, 
        {'DSPC':21, 'ffa24':43}, {'DSPC':21, 'oh12':43}, {'DSPC':21, 'oh16':43},
        {'DSPC':21, 'oh24':43}]

for composition_i in all_compositions:
    # Initialize data structures
    composition_data = {}
    composition_data['oh_tau_mean'] = []
    composition_data['hh_tau_mean'] = []
    composition_data['dip_tau_mean'] = []

    composition_data['oh_tau_std'] = []
    composition_data['hh_tau_std'] = []
    composition_data['dip_tau_std'] = []

    composition_data['unit'] = 'ps'

    for component, number in composition_i.items():
        composition_data[component] = number
        for i, name in enumerate(index.keys()):
            if index[name]['components'] == composition_i:
                os.chdir(os.path.join(curr_dir, name))
                oh_corr = np.loadtxt('oh_corr.dat')
                hh_corr = np.loadtxt('hh_corr.dat')
                dip_corr = np.loadtxt('dip_corr.dat')
                dt = oh_corr[1,1] - oh_corr[0,1]
                composition_data['oh_tau_mean'].append(np.sum(dt*oh_corr[:,2]))
                composition_data['hh_tau_mean'].append(np.sum(dt*hh_corr[:,2]))
                composition_data['dip_tau_mean'].append(np.sum(dt*dip_corr[:,2]))

    composition_data['oh_tau_std'] = np.std(composition_data['oh_tau_mean']) 
    composition_data['oh_tau_mean'] = np.mean(composition_data['oh_tau_mean']) 

    composition_data['hh_tau_std'] = np.std(composition_data['hh_tau_mean']) 
    composition_data['hh_tau_mean'] = np.mean(composition_data['hh_tau_mean']) 

    composition_data['dip_tau_std'] = np.std(composition_data['dip_tau_mean']) 
    composition_data['dip_tau_mean'] = np.mean(composition_data['dip_tau_mean']) 

    summary_df = summary_df.append(composition_data, ignore_index=True)


os.chdir(curr_dir)
summary_df.to_csv('wor.csv')
