import logging
import os
import json
import pdb
import subprocess
import itertools

import mdtraj

import gmx_ndx_functions

################
## Computing residue-water hydrogen bonds
## Not just not-water + water hbonds,
## but DSPC-water or OH18-water hbonds
##############
logging.basicConfig(filename='specific_hbonds.log',level=logging.DEBUG)

curr_dir = os.getcwd()
index = json.load(open('combined_index.json' ,'r'))

index_selections = gmx_ndx_functions.get_index_selections()

for comp in index.keys():
    logging.info(comp)
    os.chdir(os.path.join(curr_dir, comp))
    gmx_ndx_functions.write_ndx()
    components = [key for key in index[comp]['components'].keys()]
    components.append('HOH')
    for first, second in itertools.combinations_with_replacement(components,2):
        logging.info('{0}-{1}'.format(first,second))
        with open('{0}-{1}_hbond.log'.format(first, second), 'w') as f:
            p = subprocess.Popen(
                    ('echo {0} {1} | gmx hbond -f npt_80-100ns.xtc '
                    '-s npt.tpr -n hbond.ndx -life {0}-{1}_life.xvg '
                    '-num {0}-{1}_num.xvg'.format(first, second)),
                shell=True, stdout=f, stderr=f)
            p.wait()
