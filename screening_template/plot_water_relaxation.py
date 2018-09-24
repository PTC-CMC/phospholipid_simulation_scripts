import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import pdb
import plot_ay
plot_ay.setDefaults()

df = pd.read_csv('wor.csv')

bar_width = 0.5

############################
### O-H vector correlation ###
############################
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 43.0) & (df['DSPC'] == 21.0)]
    ffabar_vals.append(thing['oh_tau_mean'].values[0])
    ffayerr_vals.append(thing['oh_tau_std'].values[0])


ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 43.0) & (df['DSPC'] == 21.0)]
    ohbar_vals.append(thing['oh_tau_mean'].values[0])
    ohyerr_vals.append(thing['oh_tau_std'].values[0])


oh_x_vals = 1.5*np.arange(0, len(ohbar_vals), dtype=int) + 0.5
ffa_x_vals = 1.5*np.arange(0, len(ffabar_vals), dtype=int) + 1
xtick_vals = 1.5*np.arange(0, len(ffabar_vals), dtype=int) + 0.75

oh_bars = ax.bar(oh_x_vals, ohbar_vals, bar_width,
                    yerr=ohyerr_vals, label='67% OH, 33% DSPC')
ffa_bars = ax.bar(ffa_x_vals, ffabar_vals, bar_width,
                    yerr=ffayerr_vals, label='67% FFA, 33% DSPC')

ax.set_xticks(xtick_vals)
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("O-H relaxation (ps)")
ax.set_xlabel("Tail length")
plot_ay.tidyUp(fig, ax, legendArgs={'loc':3},
        tightLayoutArgs={})
fig.savefig("oh_relaxation.png")
plt.close(fig)

############################
### H-H vector correlation ###
############################
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 43.0) & (df['DSPC'] == 21.0)]
    ffabar_vals.append(thing['hh_tau_mean'].values[0])
    ffayerr_vals.append(thing['hh_tau_std'].values[0])


ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 43.0) & (df['DSPC'] == 21.0)]
    ohbar_vals.append(thing['hh_tau_mean'].values[0])
    ohyerr_vals.append(thing['hh_tau_std'].values[0])


oh_x_vals = 1.5*np.arange(0, len(ohbar_vals), dtype=int) + 0.5
ffa_x_vals = 1.5*np.arange(0, len(ffabar_vals), dtype=int) + 1
xtick_vals = 1.5*np.arange(0, len(ffabar_vals), dtype=int) + 0.75

oh_bars = ax.bar(oh_x_vals, ohbar_vals, bar_width,
                    yerr=ohyerr_vals, label='67% OH, 33% DSPC')
ffa_bars = ax.bar(ffa_x_vals, ffabar_vals, bar_width,
                    yerr=ffayerr_vals, label='67% FFA, 33% DSPC')

ax.set_xticks(xtick_vals)
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("H-H relaxation (ps)")
ax.set_xlabel("Tail length")
plot_ay.tidyUp(fig, ax, legendArgs={'loc':3},
        tightLayoutArgs={})
fig.savefig("hh_relaxation.png")
plt.close(fig)


############################
### Dipole vector correlation ###
############################
fig, ax = plt.subplots(1,1)
ffabar_vals = []
ffayerr_vals = []
ffa_components = ['ffa12', 'ffa16', 'ffa24']
for component in ffa_components:
    thing = df.loc[(df[component] == 43.0) & (df['DSPC'] == 21.0)]
    ffabar_vals.append(thing['dip_tau_mean'].values[0])
    ffayerr_vals.append(thing['dip_tau_std'].values[0])


ohbar_vals = []
ohyerr_vals = []
oh_components = ['oh12', 'oh16', 'oh24']
for component in oh_components:
    thing = df.loc[(df[component] == 43.0) & (df['DSPC'] == 21.0)]
    ohbar_vals.append(thing['dip_tau_mean'].values[0])
    ohyerr_vals.append(thing['dip_tau_std'].values[0])


oh_x_vals = 1.5*np.arange(0, len(ohbar_vals), dtype=int) + 0.5
ffa_x_vals = 1.5*np.arange(0, len(ffabar_vals), dtype=int) + 1
xtick_vals = 1.5*np.arange(0, len(ffabar_vals), dtype=int) + 0.75

oh_bars = ax.bar(oh_x_vals, ohbar_vals, bar_width,
                    yerr=ohyerr_vals, label='67% OH, 33% DSPC')
ffa_bars = ax.bar(ffa_x_vals, ffabar_vals, bar_width,
                    yerr=ffayerr_vals, label='67% FFA, 33% DSPC')

ax.set_xticks(xtick_vals)
ax.set_xticklabels(['12', '16', '24'])
ax.set_ylabel("Dipole relaxation (ps)")
ax.set_xlabel("Tail length")
plot_ay.tidyUp(fig, ax, legendArgs={'loc':3},
        tightLayoutArgs={})
fig.savefig("dip_relaxation.png")
plt.close(fig)

