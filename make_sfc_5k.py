import numpy as np
import subprocess as sp
import os
import datetime
from mpi4py import MPI


# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

out_dir = '/home/tsw35/soteria/clubb/data/surfaces_5k'
nate_dir = '/stor/soteria/hydro/private/nc153/projects/CLASP/LASSO/Postprocess4LES/workspace'
day_dir  = '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
day_list=[]
proj = '"+proj=aea +lat_1=34.20 +lat_2=39.20 +lat_0=36.7 +lon_0=-97.6 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"'
i=0
for file in os.listdir(nate_dir+'/sh/')[rank::size]:
    shloc=nate_dir+'/sh/'+file+'/sim.vrt'
    lhloc=nate_dir+'/lh/'+file+'/sim.vrt'
    lwloc=nate_dir+'/lwup/'+file+'/sim.vrt'
    cmd = 'gdalwarp -r average -tr 5000 5000 -te -50000 -50000 50000 50000 -t_srs '+proj+' '+shloc+' '+out_dir+'/sh/'+file+'.tif'
    #sp.run(cmd,shell=True)
    cmd = 'gdalwarp -r average -tr 5000 5000 -te -50000 -50000 50000 50000 -t_srs '+proj+' '+lhloc+' '+out_dir+'/lh/'+file+'.tif'
    #sp.run(cmd,shell=True)
    cmd = 'gdalwarp -r average -tr 5000 5000 -te -50000 -50000 50000 50000 -t_srs '+proj+' '+lwloc+' '+out_dir+'/lw/'+file+'.tif'
    sp.run(cmd,shell=True)


