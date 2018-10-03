import os
import json
import pdb
import subprocess
from collections import OrderedDict
import numpy as np

import mdtraj

import bilayer_analysis_functions

curr_dir = os.getcwd()
index = json.load(open('index.txt' ,'r'))

for comp in index.keys():
    print(comp)
    os.chdir(os.path.join(curr_dir, comp))

    traj = mdtraj.load('npt_80-100ns.xtc', top='npt.gro')
    idig, idig_avg, idig_std, = bilayer_analysis_functions.calc_interfacial_interdigitation(traj)
    np.savetxt('interface_idig.dat', np.array(idig))
    data = OrderedDict()
    data['interface_idig_mean'] = idig_avg._value
    data['interface_idig_std'] = idig_std._value
    data['interface_idig_unit'] = str(idig_avg.unit)
    with open('interface_idig.txt', 'w') as f:
        json.dump(data, f, indent=2)

