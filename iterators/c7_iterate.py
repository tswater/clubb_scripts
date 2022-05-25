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
k = 1
c7s = [0.1,0.3,0.4,0.5,0.6,0.7,0.9]
clubb_dir= '/home/tsw35/tyche/clubb/'
c7_dir = 'cpl_c7_testing/'

try:
    os.mkdir(clubb_dir+c7_dir)
except:
    pass

for l in c7s[rank::size]:
    cmd1 = 'python cpl_ant.py -v '+str(l)+' -i '+c7_dir+'v_'+str(l)
    print('#######\n'+cmd1+'\n'+'rank: '+str(rank)+'\n#######')
    cmd2 = 'python cpl_agg.py -i '+clubb_dir+c7_dir+'v_'+str(l)+'/'
    subprocess.run(cmd1,shell=True)
    print('#######\n'+cmd2+'\n'+'rank: '+str(rank)+'\n#######')
    subprocess.run(cmd2,shell=True)

