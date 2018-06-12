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
all_compositions = [{'DSPC':32, 'ffa12':32}, {'DSPC':32, 'ffa16':32}, 
        {'DSPC':32, 'ffa24':32}, {'DSPC':32, 'oh12':32}, {'DSPC':32, 'oh16':32},
        {'DSPC':32, 'oh24':32}]
for composition_i in all_compositions:
    # Initialize data structures
    composition_data = {}
    composition_data['APL_mean'] = []
    composition_data['APL_std'] = []
    composition_data['APL_unit'] = []

    composition_data['APT_mean'] = []
    composition_data['APT_std'] = []
    composition_data['APT_unit'] = []

    composition_data['H_mean'] = []
    composition_data['H_std'] = []
    composition_data['H_unit'] = []

    composition_data['Tilt_mean'] = []
    composition_data['Tilt_std'] = []
    composition_data['Tilt_unit'] = []

    composition_data['S2_mean'] = []
    composition_data['S2_std'] = []

    composition_data['Idig_mean'] = []
    composition_data['Idig_std'] = []
    composition_data['Idig_unit'] = []

    for component, number in composition_i.items():
        composition_data[component] = number
        composition_data[component + "_offset_mean"] = []
        composition_data[component + "_offset_std"] = []
        composition_data[component + "_offset_unit"] = []

    # Aggregate all the associated raw data from each simulation
    for i, name in enumerate(index.keys()):
        os.chdir(os.path.join(curr_dir, name))
        sim_data = json.load(open('data.txt','r'), object_pairs_hook=OrderedDict)

        if index[name]['components'] == composition_i:
            composition_data['APL_mean'].append(sim_data['APL']['mean'])
            composition_data['APL_std'].append(sim_data['APL']['std'])
            composition_data['APL_unit'] = sim_data['APL']['unit']

            composition_data['APT_mean'].append(sim_data['APT']['mean'])
            composition_data['APT_std'].append(sim_data['APT']['std'])
            composition_data['APT_unit'] = sim_data['APT']['unit']

            composition_data['H_mean'].append(sim_data['Bilayer Height']['mean'])
            composition_data['H_std'].append(sim_data['Bilayer Height']['std'])
            composition_data['H_unit'] = sim_data['Bilayer Height']['unit']

            composition_data['Tilt_mean'].append(sim_data['Tilt Angle']['Bilayer']['mean'])
            composition_data['Tilt_std'].append(sim_data['Tilt Angle']['Bilayer']['std'])
            composition_data['Tilt_unit'] = sim_data['Tilt Angle']['unit']

            composition_data['S2_mean'].append(sim_data['S2']['mean'])
            composition_data['S2_std'].append(sim_data['S2']['std'])

            composition_data['Idig_mean'].append(sim_data['Interdigitation']['mean'])
            composition_data['Idig_std'].append(sim_data['Interdigitation']['std'])
            composition_data['Idig_unit'] = sim_data['Interdigitation']['unit']

            for component in composition_i.keys():
                composition_data[component + "_offset_mean"].append(
                                        sim_data['Offset'][component]['mean'])
                composition_data[component + "_offset_std"].append(
                                        sim_data['Offset'][component]['std'])
                composition_data[component + "_offset_unit"] = sim_data['Offset']['unit']
    all_df = all_df.append(composition_data, ignore_index=True)

    # Now compute averages and propagate error
    ufloats = [ufloat(val,std) for val,std in zip(composition_data['APL_mean'], 
                                                  composition_data['APL_std'])]
    foo = np.mean(ufloats)
    composition_data['APL_mean'] = foo.n
    composition_data['APL_std'] = foo.s

    #composition_data['APL_mean'] = np.mean(composition_data['APL_mean'])
    #composition_data['APL_std'] = np.sqrt(sum(err**2 for err in composition_data['APL_std']))

    ufloats = [ufloat(val,std) for val,std in zip(composition_data['APT_mean'], 
                                                  composition_data['APT_std'])]
    foo = np.mean(ufloats)
    composition_data['APT_mean'] = foo.n
    composition_data['APT_std'] = foo.s
    #composition_data['APT_mean'] = np.mean(composition_data['APT_mean'])
    #composition_data['APT_std'] = np.sqrt(sum(err**2 for err in composition_data['APT_std']))

    ufloats = [ufloat(val,std) for val,std in zip(composition_data['H_mean'], 
                                                  composition_data['H_std'])]
    foo = np.mean(ufloats)
    composition_data['H_mean'] = foo.n
    composition_data['H_std'] = foo.s
    #composition_data['H_mean'] = np.mean(composition_data['H_mean'])
    #composition_data['H_std'] = np.sqrt(sum(err**2 for err in composition_data['H_std']))

    ufloats = [ufloat(val,std) for val,std in zip(composition_data['Tilt_mean'], 
                                                 composition_data['Tilt_std'])]
    foo = np.mean(ufloats)
    composition_data['Tilt_mean'] = foo.n
    composition_data['Tilt_std'] = foo.s
    #composition_data['Tilt_mean'] = np.mean(composition_data['Tilt_mean'])
    #composition_data['Tilt_std'] = np.sqrt(sum(err**2 for err in composition_data['Tilt_std']))

    ufloats = [ufloat(val,std) for val,std in zip(composition_data['S2_mean'], 
                                                  composition_data['S2_std'])]
    foo = np.mean(ufloats)
    composition_data['S2_mean'] = foo.n
    composition_data['S2_std'] = foo.s
    #composition_data['S2_mean'] = np.mean(composition_data['S2_mean'])
    #composition_data['S2_std'] = np.sqrt(sum(err**2 for err in composition_data['S2_std']))

    ufloats = [ufloat(val,std) for val,std in zip(composition_data['Idig_mean'], 
                                                  composition_data['Idig_std'])]
    foo = np.mean(ufloats)
    composition_data['Idig_mean'] = foo.n
    composition_data['Idig_std'] = foo.s
    #composition_data['Idig_mean'] = np.mean(composition_data['Idig_mean'])
    #composition_data['Idig_std'] = np.sqrt(sum(err**2 for err in composition_data['Idig_std']))

    for component, number in composition_i.items():
        ufloats = [ufloat(val,std) for val,std in 
                    zip(composition_data[component + '_offset_mean'], 
                        composition_data[component + '_offset_std'])]
        foo = np.mean(ufloats)
        composition_data[component + '_offset_mean'] = foo.n
        composition_data[component + '_offset_std'] = foo.s
        #composition_data[component + "_offset_mean"] = np.mean(composition_data[component + "_offset_mean"])
        #composition_data[component + "_offset_std"] = np.sqrt(sum(err**2 for err in composition_data[component+'_offset_std']))

    summary_df = summary_df.append(composition_data, ignore_index=True)

os.chdir(curr_dir)
summary_df.to_csv('summary.csv')
all_df.to_csv('all_sims.csv')
