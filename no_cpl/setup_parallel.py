# Sets up run folders, reads in surface, clusters as appropriate
# 

import subprocess
import numpy as np
import os
import datetime
import sklearn.cluster
import netCDF4 as nc
from mpi4py import MPI 
import shutil

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
clubb_dir = '/home/tsw35/tyche/clubb/' # directory of clubb run folders
blank_run = '/home/tsw35/soteria/clubb/clubb_scripts/run_/' # a clean run folder to copy
cbin_dir  = '/home/tsw35/soteria/software/CLUBB/bin/' # bin for clubb

# run options
dirname = 'conv_0v2' # name of main directory created
k       = 50 # (max) number of clusters 
conv    = True # if True, will setup multiple runs up to k clusters
nx      = 1000 # number of gridcells in each direction
dx      = 100 # resolution of surface grid in meters
stdate  = '2017-07-16T10:00:00.000'# start date iso format
enddate = '2017-07-17T03:00:00.000'# end date iso format
save_c  = True # True will save the clusters and surface grids
customk = True
skip_zm = False
skip_sf = False
clst_1  = True # if true, only cluster one time for whole time domain
utcoff  = -6
kcust   = [] 

# ---------------- #
# HELPER FUNCITONS #
# ---------------- #

#### READ IN SURFACE DATA ####
def read_sfc_data(var,nt,stdt,override='X'):
    bnm = 'jss'+var+'_bdy_02_' #base name
    varout_g = np.zeros((nt,nx,nx))
    varout_v = np.zeros((nt,nx*nx))
    if nt == 1:
        file=override
        with open(file,'r') as fp:
            varout_v = np.array([float(i) for i in fp.readlines()])
        varout_g=np.reshape(varout_v,(nx,nx))
        return varout_g,varout_v
    for t in range(nt):
        dt = stdt+datetime.timedelta(hours=t)
        tnm = dt.strftime('%Y-%m-%d-%H-%M')
        file = sfc_dir+bnm+tnm
        with open(file,'r') as fp:
            var_v = np.array([float(i) for i in fp.readlines()])
        varout_g[t,:,:] = np.reshape(var_v,(nx,nx))
        varout_v[t,:] = var_v[:]
    return varout_g,varout_v

#### NORMALIZE THE DATA ####
def normalize(data):
    out=(data-np.min(data))/(np.max(data)-np.min(data))
    return out


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

# adjust modelin
if rank == 0:
    lines=[]
    fp = open(blank_run+'input/arm_model.in','r')
    for line in fp:
        if line[0:3]=='day':
            dln ='day   = '+str(stdt.day)+'   ! Day of model start (1 to 31)'
            lines.append(dln+'\n')
        elif line[0:5]=='month':
            mln = 'month = '+str(stdt.month)+'   ! Month of model start (1 to 12)'
            lines.append(mln+'\n')
        elif line[0:4]=='year':
            yln = 'year  = '+str(stdt.year)+'  ! Year of model start'
            lines.append(yln+'\n')
        elif line[0:12]=='time_initial':
            iln = 'time_initial = '+str(t_init).ljust(9)+\
              '! Model start time [seconds since midnight UTC on start date]'
            lines.append(iln+'\n')
        elif line[0:10]=='time_final':
            fln = 'time_final = '+str(t_final).ljust(11)+\
              '! Model end time [seconds since midnight UTC on start date]'
            lines.append(fln+'\n')
        else:
            lines.append(line)
    fp.close()
    fp = open(blank_run+'input/arm_model.in','w')
    for line in lines:
        fp.write(line)
    fp.close()

# create file structure; copy run folders
w_dir = clubb_dir+dirname
if rank == 0:
    try:
        os.mkdir(w_dir)
    except:
        shutil.rmtree(w_dir)
        os.mkdir(w_dir)
comm.Barrier()

# figure out what clustering options we are doing
if conv:
    klist=list(range(1,k+1))
if customk:
    klist=kcust
elif not conv:
    klist=[1,k]
for i in list(range(len(klist)))[rank::size]:
    ki = klist[i]
    print(rank)
    try:
        os.mkdir(w_dir+'/k_'+str(ki))
    except OSError as error:
        print(error)
    for j in range(1,ki+1):
        print(w_dir+'/k_'+str(ki)+'/c_'+str(j))
        shutil.copytree(blank_run,w_dir+'/k_'+str(ki)+'/c_'+str(j))
        bindir = w_dir+'/k_'+str(ki)+'/c_'+str(j)+'/bin/'
        os.remove(bindir+'clubb_standalone')
        os.symlink(cbin_dir+'clubb_standalone',bindir+'clubb_standalone')

comm.Barrier()

# Read in the surface data (thanks for the weird format)
# assume data is hourly 
print('... read in surface data',flush=True)
nt     = int((endt-stdt).seconds/3600)+1
Hg,Hv   = read_sfc_data('sh',nt,stdt)
Lg,Lv   = read_sfc_data('lh',nt,stdt)
lwg,lwv  = read_sfc_data('lw',nt,stdt)
lon_g,lonv = read_sfc_data('',1,stdt,sfc_dir+'jsslongrid_02')
lat_g,latv = read_sfc_data('',1,stdt,sfc_dir+'jsslatgrid_02')
Tsfcv = (lwv/(5.67*10**(-8)))**(1/4)
Tsfcg = (lwg/(5.67*10**(-8)))**(1/4)

# ------- #
# CLUSTER #
# ------- #
# cluster based on sensible heat and then adjust the other surface var
print('CLUSTERING',flush=True)

k_masks=np.zeros((len(klist),nt,nx*nx))
for i in list(range(len(klist)))[rank::size]:
    ki = klist[i]
    print('...clustering for k = '+str(ki))
    if i > 0:
        for t in range(nt):
            if clst_1:
                continue
            Hvn   = normalize(Hv[t,:])
            latvn = normalize(latv[:])
            lonvn = normalize(lonv[:])
            X = np.reshape(Hvn,(-1,1))
            model = sklearn.cluster.KMeans(n_clusters=ki)
            # X_train, X_test = sklearn.model_selection.train_test_split(X,train_size=.1,random_state=1)
            model.fit(X)
            k_masks[i,t,:] = model.predict(X)
        if clst_1:
            Hv2 = Hv[np.mean(Hv,axis=1)>20,:]
            Hvn   = normalize(np.mean(Hv2[:,:],axis=0))
            X = np.reshape(Hvn,(-1,1))
            model = sklearn.cluster.KMeans(n_clusters=ki)
            model.fit(X)
            Y = model.predict(X)
            for t in range(nt):
                k_masks[i,t,:] = Y[:]
    else:
        k_masks[i,:,:] = np.zeros((nt,nx*nx))

comm.Barrier()

# ------------------- #
# WRITE SURFACE FILES #
# ------------------- #
# transfer clustered array to surface files
print('WRITE SURFACE FILES',flush=True)

# write initial surface files using Nate's script
for ki in klist[rank::size]:
    for j in range(1,ki+1):
        dir_old = os.getcwd()
        os.chdir(w_dir+'/k_'+str(ki)+'/c_'+str(j)+'/other_scripts')
        subprocess.run('python create_arm_data.py',shell=True)
        os.chdir(dir_old)

# overwrite arm_sfc.in
for i in list(range(len(klist)))[rank::size]:
    ki = klist[i]
    for j in range(1,ki+1):
        sfc_path = w_dir+'/k_'+str(ki)+'/c_'+str(j)+\
                   '/input/case_setups/arm_sfc.in'
        # begin overwriting file
        fp = open(sfc_path,'w')
        
        # Write the header information
        fp.write('! $Id$\n')
        fp.write('!\n')
        fp.write('! Note that T_sfc is included here,'+\
                 ' but it is not currently used. It is\n')
        fp.write('! included here to facilitate a transition '+\
                 'to using T_sfc in the future\n')
        fp.write('! if needed.\n')
        fp.write('Time[s]    latent_ht[W\m^2]   sens_ht[W\m^2]   T_sfc[K]\n')
       
        dt = 3600  # seconds

        for t in range(nt):
            ts = t_init+t*dt
            sh = np.mean(Hv[t,:][k_masks[i,t,:]==(j-1)])
            lh = np.mean(Lv[t,:][k_masks[i,t,:]==(j-1)])
            tskin = np.mean(Tsfcv[t,:][k_masks[i,t,:]==(j-1)])
            tmp = '%.2f %.2f %.2f %.2f\n' % (ts,lh,sh,tskin)
            fp.write(tmp)
        fp.close()

comm.Barrier()


# --------------- #
# OUTPUT CLUSTERS #
# --------------- #
if save_c:
    for i in list(range(len(klist)))[rank::size]:
        ki = klist[i]
        fp=nc.Dataset(w_dir+'/k_'+str(ki)+'/clusters.nc','w')
        fp.createDimension('t',size=nt)
        fp.createDimension('lat',size=nx)
        fp.createDimension('lon',size=nx)
        fp.createVariable('cluster','d',dimensions=('t','lon','lat'))
        fp.createVariable('tskin','d',dimensions=('t','lon','lat'))
        fp.createVariable('H','d',dimensions=('t','lon','lat'))
        fp.createVariable('LE','d',dimensions=('t','lon','lat'))
        fp['tskin'][:]=Tsfcg[:]
        fp['cluster'][:]=np.reshape(k_masks[i,:,:],(nt,nx,nx))
        fp['H'][:]=Hg[:]
        fp['LE'][:]=Lg[:]
        fp.close()


# ---------------- #
# RUN the MODEL(S) #
# ---------------- #
print('RUNNING MODEL',flush=True)
for i in list(range(len(klist)))[rank::size]:
    ki = klist[i]
    for j in range(1,ki+1):
        rfile = w_dir+'/k_'+str(ki)+'/c_'+str(j)+'/run_scripts/run_scm.bash'
        subprocess.run('./'+rfile+' arm',shell=True)

comm.Barrier()

# ---------------- #
# AGGREGATE OUTPUT #
# ---------------- #
print('AGGREGATE',flush=True)
btsmp = [] # timestamp for surface info
for t in range(nt):
    btsmp.append(t_init+t*3600)
times = nc.Dataset(w_dir+'/k_1/c_1/output/arm_zm.nc','r')['time'][:]
alts  = nc.Dataset(w_dir+'/k_1/c_1/output/arm_zm.nc','r')['altitude'][:]

keep_var = ['thlp2_ta','thlp2_tp','thlp2_dp1','thlp2_dp2','thlp2_ma','thlp2_forcing','thlp2']

for i in list(range(len(klist)))[rank::size]:
    ki = klist[i]

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
                if ki == 1:
                    print(str(var)+'  '+str(np.mean(wts2[:,j-1])))
                if var in ['altitude','time','longitude','latitude']:
                    continue
                fp_out[var][:,0,0,0]=fp_out[var][:,0,0,0]+\
                                     fp_in[var][:,0,0,0]*wts2[:,j-1]
                if ki == 1:
                    print(str(np.mean(fp_out[var][:]))+'   '+\
                          str(np.mean(fp_in[var][:,0,0,0]*wts2[:,j-1])))

