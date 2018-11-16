import os
import pdb
import json
import numpy as np
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd

#####################
### process trials of lipid rotations 
##############

index = json.load(open('newindex.txt', 'r'))
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()


all_compositions = [{'DSPC':21, 'ffa12':43}, {'DSPC':21, 'ffa16':43}, 
        {'DSPC':21, 'ffa24':43}, {'DSPC':21, 'oh12':43}, {'DSPC':21, 'oh16':43},
        {'DSPC':21, 'oh24':43}]


for composition_i in all_compositions:
    # Initialize data structures
    composition_data = {}
    composition_data['lip_rot_mean'] = []
    composition_data['lip_rot_std'] = []

    for component, number in composition_i.items():
        composition_data[component] = number
        for i, name in enumerate(index.keys()):
            if index[name]['components'] == composition_i:
                os.chdir(os.path.join(curr_dir, name))
                rotations = json.load(open('rotations.json', 'r'))
                composition_data['lip_rot_mean'].append(rotations['mean'])
                composition_data['lip_rot_std'].append(rotations['std'])


    ufloats = [ufloat(val,std) for val,std in zip(composition_data['lip_rot_mean'], 
                                                  composition_data['lip_rot_std'])]
    foo = np.mean(ufloats)
    composition_data['lip_rot_mean'] = foo.n
    composition_data['lip_rot_std'] = foo.s

    summary_df = summary_df.append(composition_data, ignore_index=True)


os.chdir(curr_dir)
summary_df.to_csv('lipid_rotation.csv')
