import os
import itertools
import glob
import pdb
import json
import numpy as np
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd
index = json.load(open('combined_index.json', 'r'))

curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()

compositions = [['DSPC', 'oh18'], ['DSPC', 'oh20'], ['DSPC', 'oh22']]


for pair in compositions:
    # Find simulations associated with these components
    respective_sims = [folder for folder in index.keys() if set(index[folder]['components'].keys()) == set(pair)]

    # Define all hbonding groups
    hbond_groups = [*pair, 'HOH']
    all_hbond_pairs = [*itertools.combinations_with_replacement(hbond_groups,2)]
    all_hbond_pairs.remove(('DSPC', 'DSPC'))

    # Initialize empty dictionaries to hold info from all sims
    composition_data = {}
    
    for first,second in all_hbond_pairs:
        composition_data['{0}-{1}_hbnum_mean'.format(first,second)] = []
        composition_data['{0}-{1}_hbnum_std'.format(first,second)] = []
        composition_data['{0}-{1}_hblife_mean'.format(first,second)] = []
        composition_data['{0}-{1}_hblife_std'.format(first,second)] = []
        composition_data['{0}-{1}_hblife_unit'.format(first,second)] = 'ps'

    for sim in respective_sims:
        for c in pair:
            composition_data[c] = index[sim]['components'][c]

        os.chdir(os.path.join(curr_dir, sim))
        for first,second in all_hbond_pairs:
            if os.path.isfile('{0}-{1}_life.xvg'.format(first, second)):
                life_file = np.loadtxt('{0}-{1}_life.xvg'.format(first, second),
                        comments=['@', '#'])
                num_file = np.loadtxt('{0}-{1}_num.xvg'.format(first, second),
                        comments=['@', '#'])
                calc_stuff=True
            elif os.path.isfile('{1}-{0}_life.xvg'.format(first, second)):
                life_file = np.loadtxt('{1}-{0}_life.xvg'.format(first, second),
                        comments=['@', '#'])
                num_file = np.loadtxt('{1}-{0}_num.xvg'.format(first, second),
                        comments=['@', '#'])
                calc_stuff=True
            else:
                calc_stuff=False

            if calc_stuff:
                composition_data['{0}-{1}_hbnum_mean'.format(first,second)].append(
                        np.mean(num_file[:,1]))
                composition_data['{0}-{1}_hbnum_std'.format(first,second)].append(
                        np.std(num_file[:,1]))
                composition_data['{0}-{1}_hblife_mean'.format(first,second)].append(
                        10*np.sum(life_file[:,2]))

    # After gathering all raw, sim data, propagate error
    for first, second in all_hbond_pairs:
        ufloats = [ufloat(val,std) for val,std 
            in zip(composition_data['{0}-{1}_hbnum_mean'.format(first, second)], 
                composition_data['{0}-{1}_hbnum_std'.format(first, second)])]
        foo = np.mean(ufloats)
        composition_data['{0}-{1}_hbnum_mean'.format(first, second)] = foo.n
        composition_data['{0}-{1}_hbnum_std'.format(first, second)] = foo.s

        composition_data['{0}-{1}_hblife_std'.format(first, second)] = np.std(
                composition_data['{0}-{1}_hblife_mean'.format(first, second)]) / len(
                    composition_data['{0}-{1}_hblife_mean'.format(first, second)])
        composition_data['{0}-{1}_hblife_mean'.format(first, second)] = np.mean(
                composition_data['{0}-{1}_hblife_mean'.format(first, second)])

    summary_df = summary_df.append(composition_data, ignore_index=True)
os.chdir(curr_dir)
summary_df.to_csv('specific_hbonds.csv')

