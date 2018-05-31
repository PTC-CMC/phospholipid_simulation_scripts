import os
import subprocess
import json

####################
## Remove old gmx files
##################

curr_dir = os.getcwd()
all_dirs = [thing for thing in os.listdir() if os.path.isdir(thing)]
index = json.load(open('index.txt','r'))
for folder in index.keys():
    os.chdir(os.path.join(curr_dir, folder))
    p = subprocess.Popen('rm "#"*', shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
    p.wait()
