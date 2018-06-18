import os
import json
import numpy as np
import operations
import script_utils

################
## Script to iterate through folders, write production files, and submit
################

index = json.load(open('index.txt','r'))
n_nodes = 3
jobid = [None] * n_nodes
curr_dir = os.getcwd()
for i, name in enumerate(index.keys()):
    os.chdir(os.path.join(curr_dir, name))
    lines = operations.write_production_lines(filename='npt')
    with open('production.pbs', 'w') as f:
        body = 'cd {}\n'.format(os.getcwd())
        body += lines
        script_utils.write_rahman_script(f, jobname='{}_production'.format(name), body=body)
    jobid = operations.submit_job('production.pbs', jobid, n_nodes, i)
