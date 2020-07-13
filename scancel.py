#!/usr/bin/env python
'''
python scripts to cancel a series of jobs on Slurm system
Zhihao Cui zcui@caltech.edu
'''

import sys
import subprocess

id1 = int(sys.argv[1])
id2 = int(sys.argv[2])

for i in range(id1, id2 + 1):
    subprocess.call('scancel %d'%(i), shell = True)
