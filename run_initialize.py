import numpy as np
import os
import subprocess
import json
from collections import OrderedDict

import mbuild as mb

import scripts.bilayer as bilayer

# Import statements for molecule prototypes
import atomistic.dspc.DSPC as DSPC
import atomistic.c24oh.oh24 as oh24
import atomistic.c16oh.oh16 as oh16
import atomistic.c12oh.oh12 as oh12
import atomistic.c24ffa.ffa24 as ffa24
import atomistic.c16ffa.ffa16 as ffa16
import atomistic.c12ffa.ffa12 as ffa12
import atomistic.tip3p.SOL as SOL

import operations
import script_utils

n_nodes = 2
jobid = [None] * n_nodes

curr_dir = os.getcwd()
path_to_ff = "#include \"/raid6/homes/ahy3nz/Programs/McCabeGroup/atomistic/forcefield.itp\""
solvent_density = 900
n_x = 8
n_y = 8
n_solvent_per_lipid = 20

if os.path.isfile('index.txt'):
    index = json.load(open('index.txt', 'r'))
else:
    index = OrderedDict()

table_of_contents = []

# C12oh
table_of_contents.append( ['11_dspc_c12oh_5-22-18a',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh12.oh12(), 32, -0.4) ]
                   ])

table_of_contents.append(['11_dspc_c12oh_5-22-18b',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh12.oh12(), 32, -0.4) ]
                   ])

table_of_contents.append(['11_dspc_c12oh_5-22-18c',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh12.oh12(), 32, -0.4) ]
                   ])

# C16oh
table_of_contents.append( ['11_dspc_c16oh_5-22-18a',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh16.oh16(), 32, -0.2) ]
                   ] )

table_of_contents.append(['11_dspc_c16oh_5-22-18b',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh16.oh16(), 32, -0.2) ]
                   ])

table_of_contents.append(['11_dspc_c16oh_5-22-18c',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh16.oh16(), 32, -0.2) ]
                   ])
# C24oh
table_of_contents.append(['11_dspc_c24oh_5-22-18a',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh24.oh24(), 32, -0.1) ]
                   ])

table_of_contents.append(['11_dspc_c24oh_5-22-18b',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh24.oh24(), 32, -0.1) ]
                   ])

table_of_contents.append(['11_dspc_c24oh_5-22-18c',
                    [ (DSPC.DSPC(), 32, 0),
                      (oh24.oh24(), 32, -0.1) ]
                   ])

# C12FFA
table_of_contents.append( ['11_dspc_c12ffa_5-22-18a',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa12.ffa12(), 32, -0.4) ]
                   ])

table_of_contents.append(['11_dspc_c12ffa_5-22-18b',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa12.ffa12(), 32, -0.4) ]
                   ])

table_of_contents.append(['11_dspc_c12ffa_5-22-18c',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa12.ffa12(), 32, -0.4) ]
                   ])

# C16FFA
table_of_contents.append( ['11_dspc_c16ffa_5-22-18a',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa16.ffa16(), 32, -0.2) ]
                   ] )

table_of_contents.append(['11_dspc_c16ffa_5-22-18b',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa16.ffa16(), 32, -0.2) ]
                   ])

table_of_contents.append(['11_dspc_c16ffa_5-22-18c',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa16.ffa16(), 32, -0.2) ]
                   ])
# C24FFA
table_of_contents.append(['11_dspc_c24ffa_5-22-18a',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa24.ffa24(), 32, -0.1) ]
                   ])

table_of_contents.append(['11_dspc_c24ffa_5-22-18b',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa24.ffa24(), 32, -0.1) ]
                   ])

table_of_contents.append(['11_dspc_c24ffa_5-22-18c',
                    [ (DSPC.DSPC(), 32, 0),
                      (ffa24.ffa24(), 32, -0.1) ]
                   ])


for i, composition in enumerate(table_of_contents):
    os.chdir(curr_dir)
    name = composition[0]
    lipid_info = composition[1]

    p = subprocess.Popen('mkdir -p {}'.format(name), shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    os.chdir(os.path.join(curr_dir, name))
    index[name] = OrderedDict()
    index[name]['components'] = OrderedDict()
    for species in lipid_info:
        index[name]['components'][species[0].name] = species[1]

    

    # Random apl [27,35] and tilt [0,20] and spin [0,10]
    apl = 0.30 + (0.08 * np.random.random())
    tilt_angle = np.deg2rad(20*np.random.random())
    random_spin = np.deg2rad(10*np.random.random())
    index[name]['initial_parameters'] = OrderedDict()
    index[name]['initial_parameters']['APL (nm^2)'] = apl
    index[name]['initial_parameters']['Tilt (deg)'] = np.rad2deg(tilt_angle)
    index[name]['initial_parameters']['Random spin (deg)'] = np.rad2deg(random_spin)
    with open('initial_parameters.dat', 'w') as f:
        f.write("APL (nm^2): {}\n".format(apl))
        f.write("Tilt (deg): {}\n".format(np.rad2deg(tilt_angle)))
        f.write("Random spin (deg): {}\n".format(np.rad2deg(random_spin)))

    system = bilayer.Bilayer(leaflet_info=lipid_info, n_x=n_x, n_y=n_y, apl=apl, 
            tilt_angle=tilt_angle, random_spin=random_spin,
            solvent=SOL.SOL(), solvent_density=solvent_density, 
            n_solvent_per_lipid=n_solvent_per_lipid)
    
    system = bilayer.translate_to_positive_octant(system)
    
    system.save('compound.gro', box=system.boundingbox, 
            overwrite=True, residues=set([p.parent.name for p in system.particles()]))
    bilayer.write_gmx_topology(system, 'compound.top', header=path_to_ff)

    # EQ sims
    p = subprocess.Popen('cp ~/Programs/setup/Bilayer/mdp_charmm/*.mdp .', shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    with open('eq.pbs', 'w') as f:
        body = 'cd {}\n'.format(os.getcwd())
        body += 'module load gromacs/5.1.4\n'
        body += operations.write_eq_lines(gro='compound.gro', top='compound.top')
        script_utils.write_rahman_script(f, jobname="{}_setup".format(name), body=body)

    jobid = operations.submit_job('eq.pbs', jobid, n_nodes, i)  


os.chdir(curr_dir)
with open('index.txt', 'w') as f:
    json.dump(index, f, indent=2)
