import numpy
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
df = pd.read_csv('summary.csv')

matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['axes.labelsize'] = 24
matplotlib.rcParams['legend.fontsize'] = 18

###########
### APL ###
###########
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing['APL_mean'].values)
    ffayerr_vals.append(thing['APL_std'].values)
    unit = unit_to_str[thing['APL_unit'].values[0]]

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing['APL_mean'].values)
    ohyerr_vals.append(thing['APL_std'].values)
    unit = unit_to_str[thing['APL_unit'].values[0]]


bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylim([28,36])
ax.set_ylabel("APL ({})".format(unit))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("APL.png")
plt.close(fig)

###########
### APT ###
###########
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing['APT_mean'].values)
    ffayerr_vals.append(thing['APT_std'].values)
    unit = unit_to_str[thing['APT_unit'].values[0]]

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing['APT_mean'].values)
    ohyerr_vals.append(thing['APT_std'].values)
    unit = unit_to_str[thing['APT_unit'].values[0]]


bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("APT ({})".format(unit))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("APT.png")
plt.close(fig)


##############
### Height ###
##############
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing['H_mean'].values)
    ffayerr_vals.append(thing['H_std'].values)
    unit = unit_to_str[thing['H_unit'].values[0]]

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing['H_mean'].values)
    ohyerr_vals.append(thing['H_std'].values)
    unit = unit_to_str[thing['H_unit'].values[0]]


bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("Height ({})".format(unit))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("Height.png")
plt.close(fig)

############
### Idig ###
############
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing['Idig_mean'].values)
    ffayerr_vals.append(thing['Idig_std'].values)
    unit = unit_to_str[thing['Idig_unit'].values[0]]

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing['Idig_mean'].values)
    ohyerr_vals.append(thing['Idig_std'].values)
    unit = unit_to_str[thing['Idig_unit'].values[0]]


bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("Interdigitation ({})".format(unit))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("Idig.png")
plt.close(fig)

############
### TIlt ###
############
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing['Tilt_mean'].values)
    ffayerr_vals.append(thing['Tilt_std'].values)
    unit = unit_to_str[thing['Tilt_unit'].values[0]]

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing['Tilt_mean'].values)
    ohyerr_vals.append(thing['Tilt_std'].values)
    unit = unit_to_str[thing['Tilt_unit'].values[0]]


bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("Tilt Angle ({})".format(unit))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("Tilt.png")
plt.close(fig)



#############
### Offst ###
#############
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing[component+'_offset_mean'].values)
    ffayerr_vals.append(thing[component+'_offset_std'].values)
    unit = unit_to_str[thing[component+'_offset_unit'].values[0]]

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing[component+'_offset_mean'].values)
    ohyerr_vals.append(thing[component+'_offset_std'].values)
    unit = unit_to_str[thing[component+'_offset_unit'].values[0]]


bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("Component Offset ({})".format(unit))
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("Offset.png")
plt.close(fig)


##########
### S2 ###
##########
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ffabar_vals.append(thing['S2_mean'].values)
    ffayerr_vals.append(thing['S2_std'].values)

ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 32.0) & (df['DSPC'] == 32.0)]
    ohbar_vals.append(thing['S2_mean'].values)
    ohyerr_vals.append(thing['S2_std'].values)


bar = ax.errorbar(range(len(ffabar_vals)), ffabar_vals, yerr=ffayerr_vals, label='FFA')
bar = ax.errorbar(range(len(ohbar_vals)), ohbar_vals, yerr=ohyerr_vals, label='OH')
coords = bar[0].get_xydata()
ax.set_xticks([coord[0] for coord in coords])
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("S2")
ax.set_xlabel("Tail length")
ax.legend()
fig.tight_layout()
fig.savefig("S2.png")
plt.close(fig)
