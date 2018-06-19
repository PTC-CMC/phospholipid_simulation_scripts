import os
import subprocess
import json

##############
## Check if files exist or not in each directory
###################

curr_dir = os.getcwd()
all_dirs = [thing for thing in os.listdir() if os.path.isdir(thing)]
index = json.load(open('index.txt','r'))
for folder in index.keys():
    os.chdir(os.path.join(curr_dir, folder))
    if not os.path.isfile('npt_500ps.tpr'):
        print(folder)
