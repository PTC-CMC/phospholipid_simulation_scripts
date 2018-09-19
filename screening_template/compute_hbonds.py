import os
import json
import pdb
import subprocess

import mdtraj

import gmx_ndx_functions

curr_dir = os.getcwd()
index = json.load(open('index.txt' ,'r'))

index_selections = gmx_ndx_functions.get_index_selections()

group_pairs = { 
        'solute-solute': (index_selections['non-water'], index_selections['non-water']),
        'solute-solvent': (index_selections['non-water'], index_selections['HOH']),
        'solvent-solvent': (index_selections['HOH'], index_selections['HOH'])
            }

for comp in index.keys():
    print(comp)
    os.chdir(os.path.join(curr_dir, comp))
    gmx_ndx_functions.write_ndx()
    traj = mdtraj.load('npt.gro')
    valid_residues = set(list(a.name for a in traj.topology.residues 
            if 'HOH' not in a.name and 'SOL' not in a.name))
    with open('hbond2.log', 'w') as f:
        for key in group_pairs.keys():
            if not os.path.isfile('{}_life.xvg'.format(key)):
                p = subprocess.Popen('echo {0} {1} | gmx hbond -f npt_80-100ns.xtc '.format(group_pairs[key][0], group_pairs[key][1]) +
                        '-s npt.tpr -n hbond.ndx -life {0}_life.xvg -num {0}_num.xvg'.format(
                            key),
                    shell=True, stdout=f, stderr=f)
                p.wait()
