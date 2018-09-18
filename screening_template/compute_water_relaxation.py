import os
import json

import interfacial_water_relaxation
curr_dir = os.getcwd()

########
## using mdanalysis to look at relaxation of the various water vectors
#######

index = json.load(open('index.txt', 'r'))
for folder in index.keys():
    print(folder)
    os.chdir(os.path.join(curr_dir, folder))
    try:
        interfacial_water_relaxation.main()
    except:
        print("{} failed")
