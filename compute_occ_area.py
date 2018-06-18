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
import occupied_area

""" Iterate through all Data and Sim directories, 
assessing fractional occupied area profiles """

def main():
    index = json.load(open('index.txt' ,'r'))
    curr_dir = os.getcwd()
    for key in index.keys():
        os.chdir(os.path.join(curr_dir, key))
        if os.path.isfile('npt_80-100ns.xtc') and os.path.isfile('npt.gro'):
            traj = mdtraj.load('npt_80-100ns.xtc', top='npt.gro')
            occupied_area.compute_occupied_profile_all(traj)
    os.chdir(curr_dir)


if __name__ == "__main__":
    main()
        

