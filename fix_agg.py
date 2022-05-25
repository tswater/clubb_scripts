import os
import subprocess
from mpi4py import MPI
import numpy as np

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

clubb_dir='/home/tsw35/tyche/clubb/'

aggdirs=[]

print(str(rank)+' STARTING',flush=True)

for folder in os.listdir(clubb_dir+'circ_tune/'):
    for subdir in os.listdir(clubb_dir+'circ_tune/'+folder):
        aggdir=clubb_dir+'circ_tune/'+folder+'/'+subdir+'/'
        aggdirs.append(aggdir)
i=0
for aggdir in aggdirs[rank::size]:    
    cmd2 = 'python cpl_agg.py -i '+aggdir
    print('RANK: '+str(rank)+'  '+str(i)+'/'+str(len(aggdirs)/size)+'\n'+cmd2,flush=True)
    subprocess.run(cmd2,shell=True)
    i=i+1


