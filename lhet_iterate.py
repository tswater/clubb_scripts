import netCDF4 as nc
import shutil
import subprocess
import os
import datetime
from mpi4py import MPI

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ----------- #
# USER INPUTS #
# ----------- #
k = 2
lhets = [0,10,100,1000,10000,20000,30000,40000,50000,75000,100000]
clubb_dir= '/home/tsw35/tyche/clubb/'
lhet_dir = 'cpl_lhet_testing2/'

try:
    os.mkdir(clubb_dir+lhet_dir)
except:
    pass

for l in lhets[rank::size]:
    cmd1 = 'python cpl_ant.py -l '+str(l)+' -i '+lhet_dir+'l_'+str(l)
    print('#######\n'+cmd1+'\n'+'rank: '+str(rank)+'\n#######')
    cmd2 = 'python cpl_agg.py -i '+clubb_dir+lhet_dir+'l_'+str(l)+'/'
    subprocess.run(cmd1,shell=True)
    print('#######\n'+cmd2+'\n'+'rank: '+str(rank)+'\n#######')
    subprocess.run(cmd2,shell=True)

