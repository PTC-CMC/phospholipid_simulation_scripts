import json
import numpy as np
import mdtraj

import moving_s2

if os.path.isfile('newindex.txt'):
    index = json.load(open('newindex.txt', 'r'))
else:
    index = json.load(open('index.txt', 'r'))
for composition in index.keys():
    print(composition)
    os.chdir(os.path.join(curr_dir, ratio, composition))
    traj = mdtraj.load('npt_80-100ns.xtc', top='npt.gro')
    data = moving_s2.moving_s2_routine(traj, 
            forcefield='charmm36', window_size=3)
    # Clean up data because np arrays aren't json-serializable
    data = {key1:{key2: val2.tolist() 
        for (key2, val2) in val1.items()} 
        for (key1, val1) in data.items()}
    with open('moving_s2.json', 'w') as f:
        json.dump(data, f, indent=2)
