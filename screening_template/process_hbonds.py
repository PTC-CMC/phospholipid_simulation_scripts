import os
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
    # Initialize data structures
    composition_data = {}
    composition_data['solute-solute_hbnum_mean'] = []
    composition_data['solute-solute_hbnum_std'] = []
    composition_data['solute-solute_hblife_mean'] = []
    composition_data['solute-solute_hblife_std'] = []
    composition_data['solute-solute_hblife_unit'] = 'ps'

    composition_data['solute-solvent_hbnum_mean'] = []
    composition_data['solute-solvent_hbnum_std'] = []
    composition_data['solute-solvent_hblife_mean'] = []
    composition_data['solute-solvent_hblife_std'] = []
    composition_data['solute-solvent_hblife_unit'] = 'ps'

    composition_data['solvent-solvent_hbnum_mean'] = []
    composition_data['solvent-solvent_hbnum_std'] = []
    composition_data['solvent-solvent_hblife_mean'] = []
    composition_data['solvent-solvent_hblife_std'] = []
    composition_data['solvent-solvent_hblife_unit'] = 'ps'

    for component, number in composition_i.items():
        composition_data[component] = number

    # Aggregate all the associated raw data from each simulation
    for i, name in enumerate(index.keys()):
        if index[name]['components'] == composition_i:
            os.chdir(os.path.join(curr_dir, name))
            solu_solu_hbnum = np.loadtxt('solute-solute_num.xvg', 
                    comments=['@', '#'])
            solu_solu_life = np.loadtxt('solute-solute_life.xvg', 
                    comments=['@', '#'])

            solu_solv_hbnum = np.loadtxt('solute-solvent_num.xvg', 
                    comments=['@', '#'])
            solu_solv_life = np.loadtxt('solute-solvent_life.xvg', 
                    comments=['@', '#'])

            solv_solv_hbnum = np.loadtxt('solvent-solvent_num.xvg', 
                    comments=['@', '#'])
            solv_solv_life = np.loadtxt('solvent-solvent_life.xvg', 
                    comments=['@', '#'])


            composition_data['solute-solute_hbnum_mean'].append(
                    np.mean(solu_solu_hbnum[:,1]))
            composition_data['solute-solute_hbnum_std'].append(
                    np.std(solu_solu_hbnum[:,1]))
            composition_data['solute-solute_hblife_mean'].append(
                    10*np.sum(solu_solu_life[:,2]))

            composition_data['solute-solvent_hbnum_mean'].append(
                    np.mean(solu_solv_hbnum[:,1]))
            composition_data['solute-solvent_hbnum_std'].append(
                    np.std(solu_solv_hbnum[:,1]))
            composition_data['solute-solvent_hblife_mean'].append(
                    10*np.sum(solu_solv_life[:,2]))

            composition_data['solvent-solvent_hbnum_mean'].append(
                    np.mean(solv_solv_hbnum[:,1]))
            composition_data['solvent-solvent_hbnum_std'].append(
                    np.std(solv_solv_hbnum[:,1]))
            composition_data['solvent-solvent_hblife_mean'].append(
                    10*np.sum(solv_solv_life[:,2]))




    #all_df = all_df.append(composition_data, ignore_index=True)

    # Now compute averages and propagate error
    ufloats = [ufloat(val,std) for val,std 
            in zip(composition_data['solute-solute_hbnum_mean'], 
                composition_data['solute-solute_hbnum_std'])]
    foo = np.mean(ufloats)
    composition_data['solute-solute_hbnum_mean'] = foo.n
    composition_data['solute-solute_hbnum_std'] = foo.s

    ufloats = [ufloat(val,std) for val,std 
            in zip(composition_data['solute-solvent_hbnum_mean'], 
                composition_data['solute-solvent_hbnum_std'])]
    foo = np.mean(ufloats)
    composition_data['solute-solvent_hbnum_mean'] = foo.n
    composition_data['solute-solvent_hbnum_std'] = foo.s

    ufloats = [ufloat(val,std) for val,std 
            in zip(composition_data['solvent-solvent_hbnum_mean'], 
                composition_data['solvent-solvent_hbnum_std'])]
    foo = np.mean(ufloats)
    composition_data['solvent-solvent_hbnum_mean'] = foo.n
    composition_data['solvent-solvent_hbnum_std'] = foo.s

    composition_data['solute-solute_hblife_std'] = np.std(
            composition_data['solute-solute_hblife_mean']) / len(
                    composition_data['solute-solute_hblife_mean'])
    composition_data['solute-solute_hblife_mean'] = np.mean(
            composition_data['solute-solute_hblife_mean'])

    composition_data['solute-solvent_hblife_std'] = np.std(
            composition_data['solute-solvent_hblife_mean']) / len(
                composition_data['solute-solvent_hblife_mean'])
    composition_data['solute-solvent_hblife_mean'] = np.mean(
            composition_data['solute-solvent_hblife_mean'])

    composition_data['solvent-solvent_hblife_std'] = np.std(
            composition_data['solvent-solvent_hblife_mean']) / len(
                    composition_data['solvent-solvent_hblife_mean'])
    composition_data['solvent-solvent_hblife_mean'] = np.mean(
            composition_data['solvent-solvent_hblife_mean'])

    summary_df = summary_df.append(composition_data, ignore_index=True)

os.chdir(curr_dir)
summary_df.to_csv('hbonds.csv')
#all_df.to_csv('all_sims.csv')

