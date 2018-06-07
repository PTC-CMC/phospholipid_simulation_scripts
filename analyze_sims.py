import os
import json
import numpy as np
import operations
import script_utils
import subprocess

################
## Script to iterate through folders and analyze individual simulations
################

index = json.load(open('index.txt','r'))
curr_dir = os.getcwd()
for i, name in enumerate(index.keys()):
    os.chdir(os.path.join(curr_dir, name))
    
    # Generate a gro file and pdb file at 50ns
    if not os.path.isfile('npt.gro'):
        p = subprocess.Popen('echo "0 0" | gmx trjconv -f npt.xtc -pbc mol -s npt.tpr -b 50000 -e 50000 -o npt.gro', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)  
        p.wait()
    if not os.path.isfile('npt.pdb'):
        p = subprocess.Popen('echo "0 0" | gmx trjconv -f npt.xtc -pbc mol -s npt.tpr -b 50000 -e 50000 -o npt.pdb -conect', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)  
        p.wait()

    # truncate the trajectory for last 20ns
    if not os.path.isfile('npt_30-50ns.xtc'):
        p = subprocess.Popen('echo "0 0" | gmx trjconv -f npt.xtc -pbc mol -s npt.tpr -b 30000 -e 50000 -o npt_30-50ns.xtc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)  
        p.wait()

    operations.analysis_routine('npt_30-50ns.xtc', 'npt.gro', 'npt.pdb')
