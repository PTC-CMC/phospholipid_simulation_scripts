import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import pdb

unit_to_str = { 'angstrom**2': '$\AA^2$',
                'angstrom' : '$\AA$',
                'nanometer' : 'nm',
                'nanometer**2': 'nm$^2$' ,
                'degree': 'degree'}
df = pd.read_csv('interface_idig.csv')

matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['axes.labelsize'] = 24
matplotlib.rcParams['axes.titlesize'] = 24
matplotlib.rcParams['legend.fontsize'] = 18

###########
### solute-solute hb numbers###
###########
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
oh_components = ['oh12', 'oh16','oh24']
val_matrix = np.zeros((len(ffa_components), len(oh_components)))
for i, oh in enumerate(oh_components):
    for j, ffa in enumerate(ffa_components):
        thing = df.loc[(df[oh] ==16.0) & (df['DSPC'] ==32.0) & (df[ffa] ==16.0)]
        val_matrix[i,j] = thing['interface_idig_mean'].values[0]
        unit = thing['interface_idig_unit'].values[0]

fig, ax = plt.subplots(1,1)
im = ax.imshow(val_matrix, cmap='viridis')
ax.set_xticks(range(len(ffa_components)))
ax.set_yticks(range(len(oh_components)))

#ax.set_xticklabels(ffa_components)
#ax.set_yticklabels(oh_components)
ax.set_xticklabels([thing.upper() for thing in ffa_components])
ax.set_yticklabels([thing.upper() for thing in oh_components])

ax.set_title("Interfacial Interdigitation ({}), \n50% DSPC".format(
    unit_to_str[unit]))

# Loop over data dimensions and create text annotations.
for i in range(len(oh_components)):
    for j in range(len(ffa_components)):
        text = ax.text(j, i, '{:2.2f}'.format(val_matrix[i, j]),
                           ha="center", va="center", color="w", backgroundcolor='black')

plt.colorbar(im)
fig.tight_layout()
fig.tight_layout()
fig.savefig('interfacial_idig.png')


