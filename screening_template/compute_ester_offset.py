import os
import json
import numpy as np
import operations
import script_utils
import subprocess
from collections import OrderedDict
import mdtraj
import bilayer_analysis_functions

################
## Script to iterate through folders and analyze individual simulations
################

index = json.load(open('newindex.txt','r'))
curr_dir = os.getcwd()
for i, name in enumerate(index.keys()):
    os.chdir(os.path.join(curr_dir, name))
    print(name)
    if os.path.isfile('npt_80-100ns.xtc') and os.path.isfile('npt.gro'):
        traj = mdtraj.load('npt_80-100ns.xtc', top='npt.gro')
        off_avg, off_std, off_list = bilayer_analysis_functions.calc_ester_offset(traj, blocked=True)

        np.savetxt('ester_offset.dat', np.array(off_list))

        data = OrderedDict()
        data['ester_offset_mean'] = float(off_avg._value)
        data['ester_offset_std'] = float(off_std._value)
        data['ester_offset_unit'] = str(off_avg.unit)
        with open('ester_offset.txt', 'w') as f:
            json.dump(data, f, indent=2)
    else:
        print('files not found')

