import os
import pdb
import glob
import subprocess
from multiprocessing import Pool

import math
import numpy as np
import json
import pandas as pd
import itertools
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import mdtraj
import simtk.unit as u
import bilayer_analysis_functions
import grid_analysis

""" Iterate through all Data and Sim directories, 
assessing interfacial variation using gridding """

def main():
    index = json.load(open('index.txt' ,'r'))
    curr_dir = os.getcwd()
    with Pool(16) as p:
            msr_pool = p.starmap(rough_routine, zip(itertools.repeat(curr_dir),index.keys()))

    df = pd.DataFrame(data=msr_pool, index=None, columns=['name', 'MSR_mean', 'MSR_std'])
    os.chdir(curr_dir)
    df.to_csv('roughness.csv')


def rough_routine(root_dir, sim_folder):
    """ In a simulation folder, identify the mdtraj Trajectory """
    os.chdir(os.path.join(root_dir, sim_folder))

    # Identify the proper files
    if os.path.isfile('npt.gro') and os.path.isfile('npt_80-100ns.xtc'):
        traj = mdtraj.load('npt_80-100ns.xtc', top='npt.gro')
        msr_results = bilayer_analysis_functions.analyze_simulation_interface(traj)
        msr_results['name'] = sim_folder
        return msr_results

    else:
        return None


if __name__ == "__main__":
    main()
        

