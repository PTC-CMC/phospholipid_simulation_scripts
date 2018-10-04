import os
import json
import pdb
import numpy as np

import mdtraj


curr_dir = os.getcwd()
index = json.load(open('index.txt' ,'r'))

for comp in index.keys():
    print(comp)
    os.chdir(os.path.join(curr_dir, comp))

    all_rdfs = []
    for chunk in mdtraj.iterload('npt_80-100ns.xtc', top='npt.gro', chunk=100):
        pairs = chunk.topology.select_pairs(
                'resname DSPC and (name O22 or name O32)',
                'resname HOH and name O')
        bins, rdf = mdtraj.compute_rdf(chunk, pairs, r_range=[0,2])
        all_rdfs.append(rdf)

    np.savetxt('ester-water_rdf.txt', np.column_stack([
        bins, 
        np.mean(all_rdfs, axis=0), 
        np.std(all_rdfs,axis=0)]), 
        header='r, g(r), error')

