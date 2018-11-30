import os
import json
import numpy as np

import scd_order


curr_dir = os.getcwd()
index = json.load(open('index.txt', 'r'))
    for composition in index.keys():
        print(composition)
        os.chdir(os.path.join(curr_dir, ratio, composition))
        scd_order.main(grofile='npt.gro', tprfile='npt.tpr', 
                xtcfile='npt_80-100ns.xtc')

