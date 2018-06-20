import os
import json
import numpy as np
from collections import OrderedDict
from uncertainties import ufloat
from uncertainties.umath import *
import pandas as pd
import pdb
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import bilayer_analysis_functions
################
## Script to aggregate simulation data by composition
################

index = json.load(open('index.txt','r'), object_hook=OrderedDict)
curr_dir = os.getcwd()
summary_df = pd.DataFrame()
all_df = pd.DataFrame()
all_compositions = [OrderedDict({'DSPC':32, 'ffa12':32}), 
                    OrderedDict({'DSPC':32, 'ffa16':32}), 
                    OrderedDict({'DSPC':32, 'ffa24':32}), 
                    OrderedDict({'DSPC':32, 'oh12':32}), 
                    OrderedDict({'DSPC':32, 'oh16':32}),
                    OrderedDict({'DSPC':32, 'oh24':32})]
for composition_i in all_compositions:
    all_occ_profiles = []
    for i, name in enumerate(index.keys()):

        if index[name]['components'] == composition_i:
            all_occ_profiles.append(np.loadtxt(os.path.join(name,'occupied_area_profile.dat')))
    # Since each profile will have different Z values,
    # bin the z values to group them up
    min_bin = 1000
    max_bin = 0
    for occ_profile in all_occ_profiles:
        min_bin = min(occ_profile[0,0], min_bin)
        max_bin = max(occ_profile[-1,0], max_bin)
        bin_width = occ_profile[1, 0] - occ_profile[0, 0]
    bounds = (min_bin, max_bin+0.1)
    n_bins = int(round((bounds[1] - bounds[0]) / bin_width))
    bin_edges = np.linspace(bounds[0], bounds[1], n_bins+1)
    bin_width = bin_edges[1] - bin_edges[0]

    # Cluster the values and stds based on their z_bins
    vals = [[] for x in range(len(bin_edges))]
    stds = [[] for x in range(len(bin_edges))]
    for i, occ_profile in enumerate(all_occ_profiles):
        digits = np.digitize(occ_profile[:,0], bin_edges,right=True)
        for j, digit in enumerate(digits):
            vals[digit].append(occ_profile[j, 1])
            stds[digit].append(occ_profile[j, 2])

    # Use uncertainties to compute statistics
    avg_profile = np.zeros((len(bin_edges), 2))
    for row in range(len(bin_edges)):
        ufloats = [ufloat(val, std) for val, std in zip(vals[row], stds[row])]
        if len(ufloats) > 0:
            avg_profile[row] = np.mean(ufloats).n, np.mean(ufloats).s
        else:
            avg_profile[row] = 0, 0

    bin_centers = bin_edges[:-1] + bin_width/2
    # Save and plot
    # Based on the way the binning is done, 
    # vals, stds, and avg_profile have nothing in their last element, 
    # so we can safely discard and then the bin centers and profiles have 
    # the same dimensions
    name = ""
    for key, val in sorted(composition_i.items()):
        name += '{}-{}_'.format(key,val)
    name += 'occ_profile'

    sym_profile = np.column_stack((bilayer_analysis_functions.symmetrize(avg_profile[:-1,0])[0], 
                                bilayer_analysis_functions.symmetrize(avg_profile[:-1,1])[0]))
    #np.savetxt(name+".dat", np.column_stack((bin_centers, avg_profile[:-1])))
    np.savetxt(name+".dat", np.column_stack((bin_centers, sym_profile)))
#    fig, ax = plt.subplots(1,1)
    #ax.plot(bin_centers, avg_profile[:-1, 0])
    #ax.plot(bin_centers, sym_profile[:, 0])
    #ax.set_ylim([0,1])
    #ax.fill_between(bin_centers, sym_profile[:,0] - sym_profile[:,1],
    #                sym_profile[:,0] + sym_profile[:,1], alpha=0.4)
    ##ax.fill_between(bin_centers, avg_profile[:-1,0] - avg_profile[:-1,1],
    ##                avg_profile[:-1,0] + avg_profile[:-1,1], alpha=0.4)
    #fig.savefig(name+".png")
    #plt.close(fig)
