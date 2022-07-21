#!/usr/bin/python
# Runs a crudely coupled model; 

# IMPORT 
import sys
import argparse
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


# ------------------------- #
# USER INPUTS and CONSTANTS #
# ------------------------- #

#### CONSTANTS ####
# run options
agg_zm_var = []
n_rest     = 5
delta_t    = 60 # in seconds, for CLUBB
dT_ref     = 20 # reference temperature; larger gives less exchange
dr_ref     = .015 # reference rtm; larger gives less exchange .015

# directories
sfc_dir   = '/stor/soteria/hydro/shared/sgp_surfaces/dx0100nx1000/'
clubb_dir = '../' # directory of clubb run folders
blank_run = '../run_/' # a clean run folder to copy
cbin_dir  = '/home/tsw35/soteria/software/CLUBB/bin/' # bin for clubb

#### DEFAULTS ####
k       = 2 # number of clusters
nx      = 1000 # number of gridcells in each direction
dx      = 100  # resolution of surface grid in meters
zmax    = 10000
nz      = 251
stdate  = '2017-07-16T10:00:00.000'# start date iso format
enddate = '2017-07-17T03:00:00.000'# end date iso format
dirname = 'test_cpl'

#### PARSE ####

prs = argparse.ArgumentParser(description='Short sample app')

prs.add_argument('-k', action='store', dest='k', type=int, default=k)
prs.add_argument('-n', action='store', dest='nx', type=int, default=nx)
prs.add_argument('-d', action='store', dest='dx', type=int, default=dx)
prs.add_argument('-i', action='store', dest='dir', default=dirname)
prs.add_argument('-s', action='store', dest='start', default=stdate)
prs.add_argument('-e', action='store', dest='end',default=enddate)

args = prs.parse_args()

#### FILL IN ####

k = args.k
nx = args.nx
dx = args.dx
stdate = args.start
enddate = args.end

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

#### FIND NEAREST VALUE IN AN ARRAY ####
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

#### FIND INDEX BEFORE VALUE ####
def find_previous(array,value):
    a = array-value
    a[a>=0] = -99999
    return a.argmax()

#### NORMALIZE THE DATA ####
def normalize(data):
    out=(data-np.min(data))/(np.max(data)-np.min(data))
    return out

#### READ IN THE FORCING DATA ####
def read_forcings(path,nz_=37):
    fp = open(path,'r')
    lines=fp.readlines()
    # determine number of timesteps first 
    t = 0
    for line in lines:
        if line[0]=='!':
            continue
        lsp = line.split(' ')
        if (lsp[1]==str(nz_) or (lsp[1]==(str(nz_)+'\n'))):
            t = t+1
    print(t)
    data={}
    headbool=True
    ti = -1
    alt = 0
    for line in lines:
        if line[0]=='!':
            continue
        if headbool:
            headbool = False
            lsp = line.split(' ')
            head_0 = [x for x in lsp if x != '']
            n = len(head_0)
            head=[]
            for x in head_0:
                if '\n' in x:
                    x = x[0:-1]
                head.append(x)
                data[x]=np.zeros((t,37))
            data['time']=[]
            continue
        lsp = line.split(' ')
        if (lsp[1]==str(nz_) or (lsp[1]==(str(nz_)+'\n'))):
            data['time'].append(float(lsp[0]))
            ti=ti+1
            alt=0
            continue
        for i in range(n):
            x = head[i]
            data[x][ti,alt]=float(lsp[i])
        alt=alt+1
    data['time']=np.array(data['time'])
    return data

#### WRITE FORCING DATA ####
def write_forcings(data,frc_path):
    fp = open(frc_path,'w')
    nt_ = len(data['time'])
    nz_ = len(data['Press[Pa]'][0,:])
    fp.write('! $Id$\n')
    fp.write('! The vertical coordinate is entered in the first column as height (z[m]) or pressure (Press[Pa]).\n')
    fp.write('! Temperature can be entered as liquid water potential temperature (thlm_f[K\s]), potential temperature (thm_f[K\s]), or absolute temperature(T_f[K\s]).\n')
    fp.write('! Prescribed, time-independent vertical velocity can be entered as velocity (w[m\s]) or pressure velocity (omega[Pa\s] or omega[mb\hr])\n')
    fp.write('Press[Pa]       T_f[K\s]        rtm_f[kg\kg\s]   um_ref[m\s]   vm_ref[m\s]      um_f[m\s^2]     vm_f[m\s^2]     omega[Pa\s]          ug[m\s]         vg[m\s]\n')
    for ti_ in range(nt_):
        tmp = '%.2f %d\n' % (data['time'][ti_],nz_)
        fp.write(tmp)
        for zi_ in range(nz_):
            tmp = '%.4f %.4e %.4e %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n' %\
                  (data['Press[Pa]'][ti_,zi_],data['T_f[K\s]'][ti_,zi_],\
                   data['rtm_f[kg\kg\s]'][ti_,zi_],data['um_ref[m\s]'][ti_,zi_],\
                   data['vm_ref[m\s]'][ti_,zi_],data['um_f[m\s^2]'][ti_,zi_],\
                   data['vm_f[m\s^2]'][ti_,zi_],data['omega[Pa\s]'][ti_,zi_],\
                   data['ug[m\s]'][ti_,zi_],data["vg[m\s]"][ti_,zi_])
            fp.write(tmp)
    return

# ------------- #
# INITIAL SETUP #
# ------------- #
stdt = datetime.datetime.fromisoformat(stdate)
endt = datetime.datetime.fromisoformat(enddate)

pre_dt  = datetime.datetime(stdt.year,stdt.month,stdt.day,0,0)
t_init  = (stdt-pre_dt).seconds
t_final = (endt-pre_dt).total_seconds()
tf0     = t_init+delta_t*n_rest

# create file structure; copy run folders
w_dir = clubb_dir+dirname
m_dir = w_dir+'/k_'+str(k)

if rank == 0:
    # Create or clear out directory
    try:
        os.mkdir(w_dir)
    except:
        pass
    try:
        os.mkdir(m_dir)
    except:
        shutil.rmtree(m_dir)
        os.mkdir(m_dir)
    # copy from the main run folder
    shutil.copytree(blank_run,w_dir+'/k_'+str(k)+'/c_1')
    
    # fix the arm_model.in
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
            fln = 'time_final = '+str(tf0).ljust(11)+\
              '! Model end time [seconds since midnight UTC on start date]'
            lines.append(fln+'\n')
        else:
            lines.append(line)
    fp.close()
    fp = open(m_dir+'/c_1/input/arm_model.in','w')
    for line in lines:
        fp.write(line)
    fp.close()

    # Copy and edit original to make it complete for a single, uncoupled run
    shutil.copy(m_dir+'/c_1/input/arm_model.in',m_dir+'/original_arm_model.in')
    fp = open(m_dir+'/original_arm_model.in','r')
    lines = []
    for line in fp.readlines():
        if line[0:10]=='time_final':
            fln = 'time_final = '+str(t_final).ljust(11)+\
              '! Model end time [seconds since midnight UTC on start date]'
            lines.append(fln+'\n')
        else:
            lines.append(line)
    fp.close()
    fp = open(m_dir+'/original_arm_model.in','w')
    for line in lines:
        fp.write(line)
    fp.close()
    
comm.Barrier()

for i in list(range(2,k+1))[rank::size]:
    shutil.copytree(m_dir+'/c_1',m_dir+'/c_'+str(i))
    bindir = w_dir+'/k_'+str(k)+'/c_'+str(i)+'/bin/'
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

k_masks=np.zeros((nt,nx*nx))
for i in list(range(k))[rank::size]:
    Hv2 = Hv[np.mean(Hv,axis=1)>20,:]
    Hvn   = normalize(np.mean(Hv2[:,:],axis=0))
    X = np.reshape(Hvn,(-1,1))
    model = sklearn.cluster.KMeans(n_clusters=k)
    model.fit(X)
    Y = model.predict(X)
    for t in range(nt):
        k_masks[t,:] = Y[:]

comm.Barrier()


# ------------------- #
# WRITE SURFACE FILES #
# ------------------- #
# transfer clustered array to surface files
print('WRITE SURFACE FILES',flush=True)

# write initial surface files using Nate's script
for i in list(range(k))[rank::size]:
    dir_old = os.getcwd()
    os.chdir(w_dir+'/k_'+str(k)+'/c_'+str(i+1)+'/other_scripts')
    subprocess.run('python create_arm_data_cpl.py',shell=True)
    os.chdir(dir_old)

# overwrite arm_sfc.in
for j in list(range(1,k+1)):
    sfc_path = w_dir+'/k_'+str(k)+'/c_'+str(j)+\
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
        sh = np.mean(Hv[t,:][k_masks[t,:]==(j-1)])
        lh = np.mean(Lv[t,:][k_masks[t,:]==(j-1)])
        tskin = np.mean(Tsfcv[t,:][k_masks[t,:]==(j-1)])
        tmp = '%.2f %.2f %.2f %.2f\n' % (ts,lh,sh,tskin)
        fp.write(tmp)
    fp.close()
    
    shutil.copy(sfc_path,w_dir+'/k_'+str(k)+'/arm_sfc_original.in')

comm.Barrier()

# --------------- #
# OUTPUT CLUSTERS #
# --------------- #
if rank == 0:
    fp=nc.Dataset(w_dir+'/k_'+str(k)+'/clusters.nc','w')
    fp.createDimension('t',size=nt)
    fp.createDimension('lat',size=nx)
    fp.createDimension('lon',size=nx)
    fp.createVariable('cluster','d',dimensions=('t','lon','lat'))
    fp.createVariable('tskin','d',dimensions=('t','lon','lat'))
    fp.createVariable('H','d',dimensions=('t','lon','lat'))
    fp.createVariable('LE','d',dimensions=('t','lon','lat'))
    fp['tskin'][:]=Tsfcg[:]
    fp['cluster'][:]=np.reshape(k_masks[:,:],(nt,nx,nx))
    fp['H'][:]=Hg[:]
    fp['LE'][:]=Lg[:]
    fp.close()

clst_frac = np.zeros((k,))
for i in range(k):
    clst_frac[i]=np.sum(k_masks[:]==i)/k_masks.size

nz_forcing = 37 #number of levels in forcing #FIXME 
tlist = list(range(int(t_init),int(t_final),int(round(delta_t*n_rest))))
wtzi = np.ones((len(tlist),nz_forcing,k))


# --------- #
# CORE LOOP #
# --------- #

# FIRST TIMESTEP 
for j in list(range(1,k+1))[rank::size]:
    # run the script
    rfile = w_dir+'/k_'+str(k)+'/c_'+str(j)+'/run_scripts/run_scm.bash'
    print(rfile)
    subprocess.run('./'+rfile+' arm >'+'log'+str(j),shell=True,stdin=subprocess.DEVNULL)
    print()
    print(str(j) + ' IS COMPLETE ',flush=True)
    os.mkdir(w_dir+'/k_'+str(k)+'/c_'+str(j)+'/restart')
    os.mkdir(w_dir+'/k_'+str(k)+'/c_'+str(j)+'/output/old')
comm.Barrier()


frcs_0 = read_forcings(m_dir+'/c_'+str(1)+'/input/case_setups/arm_forcings.in',nz_forcing)
p_in = frcs_0['Press[Pa]'][0,:]
t_in = frcs_0['time'][:]

for i in range(1,len(tlist)):
    # Time Naming Conventions t[NUM]i_[X] means index of t[NUM] in file [X]
    # [X] == f means forcing file, [X] == o means output file
    t3 = tlist[i] # end of previous run, start of this run
    t4 = t3+n_rest*delta_t/2 # middle of this run
    t0i_f = find_previous(t_in,t4)
    t0 = t_in[t0i_f] # lower bound of forcing file which contains t3 
    t1 = tlist[i-1] # start of the previous run
    t2 = t1+n_rest*delta_t/2 # middle of previous run
    t3 = tlist[i] # end of previous run, start of this run
    t4 = t3+n_rest*delta_t/2 # middle of this run
    t5 = t3+n_rest*delta_t # end of this run
    t6 = t_in[t0i_f+1] # upper bound of forcing file which contains t3
    
    t2i_o = int(round(n_rest/2))

    tmps = np.zeros((k,nz_forcing))
    rtms = np.zeros((k,nz_forcing))
    ums  = np.zeros((k,nz_forcing))
    vms  = np.zeros((k,nz_forcing))
    
    tmp_m =np.zeros((nz_forcing,))
    rtm_m =np.zeros((nz_forcing,))

    t3i_o = int(round(n_rest-1))
    
    # Load in the temperature and mixing ratio
    for j in range(1,k+1):
        fp = nc.Dataset(m_dir+'/c_'+str(j)+'/output/arm_zt.nc','r')
        p_out = fp['p_in_Pa'][t3i_o,:,0,0]
        
        for l in range(nz_forcing):
            zli = find_nearest(p_out,p_in[l])
            tmps[j-1,l]=fp['T_in_K'][t3i_o,zli,0,0]
            rtms[j-1,l]=fp['rtm'][t3i_o,zli,0,0]
        tmp_m[:]=tmp_m[:]+clst_frac[j-1]*tmps[j-1,:]
        rtm_m[:]=rtm_m[:]+clst_frac[j-1]*rtms[j-1,:]
        fp.close()

    # determine weights
    wts_t = np.ones((k,nz_forcing))
    wts_r = np.ones((k,nz_forcing))
    for l in range(nz_forcing):
        wts_t[:,l]=1+(tmps[:,l]-tmp_m[l])/dT_ref
        wts_t[:,l]=wts_t[:,l]/np.sum(wts_t[:,l])
        wts_r[:,l]=1+(rtms[:,l]-rtm_m[l])/dr_ref
        wts_r[:,l]=wts_r[:,l]/np.sum(wts_r[:,l])
    
    # find the index of the time before the current time
    
    for j in list(range(1,k+1))[rank::size]:
        frc_file = m_dir+'/c_'+str(j)+'/input/case_setups/arm_forcings.in'
        frcs = read_forcings(frc_file,nz_forcing)
        
        # find forcing value at center of the next run
        Et_0  = frcs['T_f[K\s]'][t0i_f,:]
        Et_6  = frcs['T_f[K\s]'][t0i_f+1,:]
        Et_4 = Et_0+(Et_6-Et_0)/(t6-t0)*(t4-t0)
        et_4 = Et_4*wts_t[j-1,:]
        et_6 = et_4+(et_4-Et_0)/(t4-t0)*(t6-t0)
        
        # write to the forcing file
        frcs['T_f[K\s]'][t0i_f+1,:]=et_6[:]
        write_forcings(frcs,frc_file)

        # copy output files to new location
        c_dir = m_dir+'/c_'+str(j)+'/'
        shutil.copy(c_dir+'output/arm_zt.nc',c_dir+'restart/arm_zt.nc')
        shutil.copy(c_dir+'output/arm_sfc.nc',c_dir+'restart/arm_sfc.nc')
        shutil.copy(c_dir+'output/arm_zm.nc',c_dir+'restart/arm_zm.nc')

        # move all the old output files out of the way
        os.rename(c_dir+'output/arm_zt.nc',c_dir+'output/old/arm_zt_'+str(i)+'.nc')
        os.rename(c_dir+'output/arm_sfc.nc',c_dir+'output/old/arm_sfc_'+str(i)+'.nc')
        os.rename(c_dir+'output/arm_zm.nc',c_dir+'output/old/arm_zm_'+str(i)+'.nc')

        # edit the model.in file
        fp = open(c_dir+'input/arm_model.in','r')
        lines= []
        for line in fp.readlines():
            if 'l_restart' in line:
                lines.append('l_restart = .true\n')
            elif 'restart_path_case' in line:
                lines.append('restart_path_case = "restart/arm"\n')
            elif 'time_restart' in line:
                lines.append('time_restart = '+str(t3)+' \n')
            elif line[0:10]=='time_final':
                fln = 'time_final = '+str(t6).ljust(11)+\
                      '! Model end time [seconds since midnight UTC on start date]'
                lines.append(fln+'\n')
            else:
                lines.append(line)
        fp.close()
        fp = open(c_dir+'input/arm_model.in','w')
        for line in lines:
            fp.write(line)
        fp.close()
    
        # run restart
        rfile = w_dir+'/k_'+str(k)+'/c_'+str(j)+'/run_scripts/run_scm.bash'
        print('STARTING RUN: '+str(j)+' at time '+str(t3)+' with process '+str(rank),flush=True)
        subprocess.run('./'+rfile+' arm',shell=True,stdout=subprocess.DEVNULL)
        print('RUN COMPLETE: '+str(j)+' at time '+str(t3)+' with process '+str(rank),flush=True)

# ---------- #
# AGG OUTPUT #
# ---------- #
# temporally aggregate small output files
# then do traditional output aggregation


        

    
    

    

