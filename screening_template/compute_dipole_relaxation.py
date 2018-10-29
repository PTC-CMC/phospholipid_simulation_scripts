import os
import json
from collections import OrderedDict
import numpy as np
from scipy.optimize import curve_fit



def stretched_exponential(x_data, a, tau, beta):
    """ Function used to get decay times for water relaxation

    Parameters
    ----------
    x_data : times 
    a : exponential prefactor 
    tau : rotational relaxation time
    beta : exponent in the exponential term


    References
    ----------
    Sciortino, P., Gallo, P., Tartaglia, P., Chen, S-H. 1996. Supercooled water
        and the kinetic glass transition
    Castrillon, S., Giovambattista, N., Aksay, I.A., Debenedetti, P. 2009. Effect
        of Surface Polarity on the Structure and Dynamics of WAter in Nanoscale
        Confinement
        """
    return a * np.exp(-((x_data/tau)**beta))


curr_dir = os.getcwd()

index = json.load(open('index.txt', 'r'))
for folder in index.keys():
    print(folder)
    os.chdir(os.path.join(curr_dir, folder))
    data = np.loadtxt('dip_corr.dat')
    times = data[:,1]
    corr = data[:,2]

    popt, pcov = curve_fit(stretched_exponential, times, corr)
    stretched_exp_params = OrderedDict()
    stretched_exp_params['A'] = popt[0]
    stretched_exp_params['tau'] = popt[1]
    stretched_exp_params['beta'] = popt[2]
    with open('dip_corr_fit.txt', 'w') as f:
        json.dump(stretched_exp_params, f, indent=2)

    
