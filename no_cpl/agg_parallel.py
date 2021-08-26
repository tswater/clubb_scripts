# Sets up run folders, reads in surface, clusters as appropriate
# 

import subprocess
import numpy as np
import os
import datetime
import sklearn.cluster
import netCDF4 as nc
from mpi4py import MPI 

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ------------------------ #
# USER INPUT and CONSTANTS #
# ------------------------ #
# user inputs for the script, directory paths, etc.

# directories
sfc_dir   = '/stor/soteria/hydro/shared/sgp_surfaces/dx0100nx1000/'
clubb_dir = '../' # directory of clubb run folders
blank_run = '../run_/' # a clean run folder to copy
cbin_dir  = '/home/tsw35/soteria/software/CLUBB/bin/' # bin for clubb

# run options
dirname = '16_7_17_conv100' # name of main directory created
k       = 10 # (max) number of clusters 
conv    = True # if True, will setup multiple runs up to k clusters
nx      = 1000 # number of gridcells in each direction
dx      = 100 # resolution of surface grid in meters
stdate  = '2017-07-16T10:00:00.000'# start date iso format
enddate = '2017-07-17T03:00:00.000'# end date iso format
save_c  = True # True will save the clusters and surface grids
skip_zm = False # Skips aggregation of zm output
skip_sf = False # Skips aggregation of surface output
customk = True
customk_arr = []
# ------------- #
# INITIAL SETUP #
# ------------- #
# setup file structure, copy files, fill arm_model.in, load surface data
print('INITIAL SETUP')
# convert the start and end dates into values required by arm_model.in
stdt = datetime.datetime.fromisoformat(stdate)
endt = datetime.datetime.fromisoformat(enddate)

pre_dt  = datetime.datetime(stdt.year,stdt.month,stdt.day,0,0)
t_init  = (stdt-pre_dt).seconds
t_final = (endt-pre_dt).total_seconds()

# create file structure; copy run folders
w_dir = clubb_dir+dirname

# figure out what clustering options we are doing
if conv:
    klist=list(range(1,k+1))
#if customk:
    klist=customk_arr 
elif not conv:
    klist=[1,k]

nt     = int((endt-stdt).seconds/3600)+1

# read in masks
k_masks=np.zeros((len(klist),nt,nx*nx))
for i in list(range(len(klist)))[rank::size]:
    ki = klist[i]
    fp = nc.Dataset(w_dir+'/k_'+str(ki)+'/clusters.nc','r')
    k_masks[i,:,:] = np.reshape(fp['cluster'][:],(nt,nx*nx))
    fp.close()

comm.Barrier()

# ---------------- #
# AGGREGATE OUTPUT #
# ---------------- #
btsmp = [] # timestamp for surface info
for t in range(nt):
    btsmp.append(t_init+t*3600)
times = nc.Dataset(w_dir+'/k_1/c_1/output/arm_zm.nc','r')['time'][:]
alts  = nc.Dataset(w_dir+'/k_1/c_1/output/arm_zm.nc','r')['altitude'][:]

keep_var = ['thlp2_ta','thlp2_tp','thlp2_dp1','thlp2_dp2','thlp2_ma','thlp2_forcing','thlp2']

for i in list(range(len(klist)))[rank::size]:
    ki = klist[i]
    print(i)

    # Create weighting scheme between clusters for all timesteps in output 
    print('first weights')
    wts_sfc = np.zeros((nt,ki))
    for j in range(ki):
        for t in range(nt):
            wts_sfc[t,j] = np.sum(k_masks[i,t,:]==j)/k_masks.shape[2]
            if ki == 1:
                print('j: '+str(j)+'   t: '+str(t)+'   wts: '+str(wts_sfc[t,j]))
    wts = np.zeros((len(times),len(alts),ki))
    wts2 = np.zeros((len(times),ki))
    idx=0
    delta = btsmp[1]-btsmp[0]
    print('make it cover full time range')
    for t in range(len(times)):
        if (btsmp[idx]+delta/2)<t:
            idx=idx+1
        for j in range(ki):
            wts[t,:,j] = wts_sfc[idx,j]
            wts2[t,j]  = wts_sfc[idx,j]

    # iterate through output for each cluster
    if not skip_zm:
        print(ki)
        for j in range(1,ki+1):
            ofile = w_dir+'/k_'+str(ki)+'/c_'+str(j)+'/output/arm_zm.nc'
            print(j)
            fp_in = nc.Dataset(ofile,'r')
            if j == 1:
                agfile = w_dir+'/k_'+str(ki)+'/agg_outzm.nc'
                fp_out=nc.Dataset(agfile,'w')
                fp_out.createDimension('time',size=len(times))
                fp_out.createDimension('latitude',size=1)
                fp_out.createDimension('longitude',size=1)
                fp_out.createDimension('altitude',size=len(alts))
                dim = ('time','altitude','latitude','longitude')
                fp_out.createVariable('time','d',dimensions=('time'))
                fp_out.createVariable('altitude','d',dimensions=('altitude'))
                #fp_out.createVariable('longitude','d',dimensions=('longitude'))
                #fp_out.createVariable('latitude','d',dimensions=('latitude'))
                fp_out['time'][:]=times[:]
                fp_out['altitude'][:]=alts[:]
                for var in keep_var:
                    fp_out.createVariable(var,'d',dimensions=dim)
                    fp_out[var][:,:,0,0]=np.zeros((len(times),len(alts)))
            for var in keep_var:
                fp_out[var][:,:,0,0]=fp_out[var][:,:,0,0]+\
                                     fp_in[var][:,:,0,0]*wts[:,:,j-1]
        
    # now again but for surface stuff
    
    if not skip_sf:
        for j in range(1,ki+1):
            ofile = w_dir+'/k_'+str(ki)+'/c_'+str(j)+'/output/arm_sfc.nc'
            fp_in = nc.Dataset(ofile,'r')
            if j == 1:
                agfile = w_dir+'/k_'+str(ki)+'/agg_outsfc.nc'
                subprocess.run('cp '+ofile+' '+agfile,shell=True)
                fp_out=nc.Dataset(agfile,'r+')
                for var in fp_out.variables:
                    if var in ['altitude','time','longitude','latitude']:
                        continue
                    fp_out[var][:,0,0,0]=np.zeros((len(times)))[:]
            for var in fp_in.variables:
                print(var)
                if ki == 1:
                    print(str(var)+'  '+str(np.mean(wts2[:,j-1])))
                if var in ['altitude','time','longitude','latitude']:
                    continue
                fp_out[var][:,0,0,0]=fp_out[var][:,0,0,0]+\
                                     fp_in[var][:,0,0,0]*wts2[:,j-1]
                if ki == 1:
                    print(str(np.mean(fp_out[var][:]))+'   '+\
                          str(np.mean(fp_in[var][:,0,0,0]*wts2[:,j-1])))
