import subprocess as sb
import numpy as np
import os


from mpi4py import MPI

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

clubdir='/home/tsw35/tyche/clubb/'

os.chdir('/home/tsw35/soteria/clubb/clubb_scripts')

r_s = [.05,.1,.5,.75,1,2,5,.05,.1,.5,.75,1,2,5,1,1]
b_s = [1,1,1,1,1,1,1,10,10,10,10,10,10,10,1,1]
k_s = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1]
f_s = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0]
for i in range(len(r_s))[rank::size]:
    diri = 'test/test%.2f_%i_%i_%i' % (r_s[i],b_s[i],k_s[i],f_s[i])
    cmd1=('time python cpl_ant.py -r %.2f -b %i -k %i -f %i -i '+diri) % (r_s[i],b_s[i],k_s[i],f_s[i])
    cmd2='python cpl_agg.py -i '+clubdir+diri+'/'
    sb.run(cmd1,shell=True)
    sb.run(cmd2,shell=True)
