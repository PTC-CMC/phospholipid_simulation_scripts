import os
import json
import numpy as np
from collections import OrderedDict
import pdb
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
################
## Script to plot occupied area fractions 
################

matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['axes.labelsize'] = 24
matplotlib.rcParams['legend.fontsize'] = 18

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
    name = ""
    for key, val in sorted(composition_i.items()):
        name += '{}-{}_'.format(key,val)
    name += 'occ_profile'

    sym_profile = np.loadtxt(name+".dat") 
    fig, ax = plt.subplots(1,1)
    ax.plot(sym_profile[:,0], sym_profile[:, 1])
    ax.set_ylim([0,1])
    ax.fill_between(sym_profile[:,0], sym_profile[:,1] - sym_profile[:,2],
                    sym_profile[:,1] + sym_profile[:,2], alpha=0.4)
    fig.savefig(name+".png")
    plt.close(fig)
