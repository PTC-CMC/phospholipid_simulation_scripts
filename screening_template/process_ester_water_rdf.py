import os
import mdtraj
import glob
import pdb
import json
import numpy as np
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd

index = json.load(open('index.txt', 'r'))
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()
all_compositions = [{'DSPC':32, 'ffa12':32}, {'DSPC':32, 'ffa16':32}, 
        {'DSPC':32, 'ffa24':32}, {'DSPC':32, 'oh12':32}, {'DSPC':32, 'oh16':32},
        {'DSPC':32, 'oh24':32}]

for composition_i in all_compositions:
    composition_data = {}
    composition_data['rdf_peak_mean'] = []
    composition_data['rdf_peak_std'] = []

    for component, number in composition_i.items():
        composition_data[component] = number

    for i, name in enumerate(index.keys()):
        if index[name]['components'] == composition_i:
            os.chdir(os.path.join(curr_dir, name))
            rdf = np.loadtxt('ester-water_rdf.txt')
            composition_data['rdf_peak_mean'].append(np.max(rdf[:,1]))

    composition_data['rdf_peak_std'] = np.std(composition_data['rdf_peak_mean'])
    composition_data['rdf_peak_mean'] = np.mean(composition_data['rdf_peak_mean'])

    summary_df = summary_df.append(composition_data, ignore_index=True)
os.chdir(curr_dir)
summary_df.to_csv('rdf_peaks.csv')

