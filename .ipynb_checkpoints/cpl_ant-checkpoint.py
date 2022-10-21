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
import shutil
from scipy.stats import beta

# ------------------------- #
# USER INPUTS and CONSTANTS #
# ------------------------- #

#### CONSTANTS ####
# run options
n_rest     = 5 # number of timesteps between restart
delta_t    = 60 # in seconds, for CLUBB  ## DOES NOT CHANGE MODEL VALUE
dsmooth    = 5 # number of levels over which to switch sign of circ flux
T0         = 300 # reference temperature; matched to CLUBB
c_r        = .00125 # factor for moisture fluxes
c_t        = .0000625 # factor for temperature fluxes

# directories
sfc_dir   = '/stor/soteria/hydro/shared/sgp_surfaces/dx0100nx1000/'
clubb_dir = '/home/tsw35/tyche/clubb/' # directory of clubb run folders
blank_run = '/home/tsw35/soteria/clubb/clubb_scripts/run_/' # a clean run folder to copy
cbin_dir  = '/home/tsw35/soteria/software/CLUBB/bin/' # bin for clubb

#### DEFAULTS ####
k       = 1 # number of clusters
nx      = 1000 # number of gridcells in each direction
dx      = 100  # resolution of surface grid in meters
dz      = 40 # dz
zmax    = 10000
dzc     = 1000
nz      = int(np.floor(zmax/dz)+1)
stdate  = '2017-07-16T12:00:00.000'# start date iso format
enddate = '2017-07-17T03:00:00.000'# end date iso format
dirname = 'test_cpl6'
l_het   =  40000 # Lengthscale of heterogeneity in meters
# FIXME as nate for a better first guess of l_het after we give him baby space
z_r     = 1000 # height where we turn off the lower near surface flux
z_r2    = 3500
z_r3    = 4500 
c_7     = -1 # default is non constant (-1), bouancy term, .3 to .8 is listed range

#### PARSE ####

prs = argparse.ArgumentParser(description='Short sample app')

prs.add_argument('-k', action='store', dest='k', type=int, default=k)
prs.add_argument('-n', action='store', dest='nx', type=int, default=nx)
prs.add_argument('-d', action='store', dest='dx', type=int, default=dx)
prs.add_argument('-i', action='store', dest='dir', default=dirname)
prs.add_argument('-s', action='store', dest='start', default=stdate)
prs.add_argument('-e', action='store', dest='end',default=enddate)
prs.add_argument('-c', action='store', dest='sfc',default=sfc_dir)
prs.add_argument('-l', action='store', dest='l_het',type=int,default=l_het)
prs.add_argument('-z', action='store', dest='z_r',type=int,default=z_r)
prs.add_argument('-r', action='store', dest='c_r',type=float,default=c_r)
prs.add_argument('-t', action='store', dest='c_t',type=float,default=c_t)
prs.add_argument('-v', action='store', dest='c_7',type=float,default=c_7)

args = prs.parse_args()

#### FILL IN ####
k       = args.k
nx      = args.nx
dx      = args.dx
stdate  = args.start
enddate = args.end
dirname = args.dir
sfc_dir = args.sfc
l_het   = args.l_het
z_r     = args.z_r
c_r     = args.c_r
c_t     = args.c_t
c_7     = args.c_7

#### EXTRAS ####
dzi = int(np.floor(dzc/dz))


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

'''
#### COMPUTE CIRCULATION FLUX (DEPRECATED) ####
# W     : width of connection (k,k)
# zs    : height in frcs space (nz)
# T     : temperature, (k,nz)
# lam   : circulating species, (k,nz)
# zr_i  : index in forcing space of switching period
# H     : surface sensible heat flux (k)
# k     : number of columns
# T0    : reference temperature
# l_het : lengthscale of heterogeneity
# F[k1,k2,z] means flux from k1 to k2. (+) is net flux out of k1
def circ_flux_dpct(W,T,lam,c_,H,V,zr_i,zr_i2=0,zr_i3=0,dz_=dz,nz_=nz,dsm=dsmooth,k=k,T0=T0,l=l_het):
    F = np.zeros((k,k,nz_))
    denom = nz_-(zr_i+dsm)
    F_sum = np.zeros((k,k))
    F_sumout = np.zeros((k,k))
    print(zr_i)
    print(zr_i2)
    print(zr_i3)
    for k1 in range(k):
        for k2 in range(k):
            if k1==k2:
                continue
            
            sgn = (H[k2]-H[k1])/np.abs(H[k2]-H[k1])
            if sgn<0:
                k_low=k2
                k_hi =k1
            else:
                k_low=k1
                k_hi =k2

            # homogeneous 2 sector case
            if zr_i2>0:
                fracsum = 0
                for i in range(1,dsm):
                    fracsum = fracsum+i/dsm
                T1m = np.mean(T[k1,0:zr_i])
                T2m = np.mean(T[k2,0:zr_i])
                ur = np.abs(T1m-T2m)/T0*9.81**(.5)*l**(.5)
                ur2 = ur*(zr_i+fracsum)/(zr_i3-zr_i2+2*fracsum)
                
                # low level flux
                F[k1,k2,0:zr_i] = c_*W[k1,k2]*dz_*ur*\
                                  np.mean(lam[k_low,0:zr_i])/(V[k_hi])*sgn
                for i in range(1,dsm):
                    F[k1,k2,zr_i+i-1]=c_*W[k1,k2]*dz_*ur*\
                                      np.mean(lam[k_low,0:zr_i])/\
                                      (V[k_hi])*sgn*(dsm-i)/(dsm)
                
                # high level flux
                F[k1,k2,zr_i2:zr_i3]=-c_*W[k1,k2]*dz_*ur2*np.mean(lam[k_hi,zr_i2:zr_i3])/(V[k_low])*sgn
                for i in range(1,dsm):
                    F[k1,k2,zr_i2-i]=-c_*W[k1,k2]*dz_*ur2*\
                                     np.mean(lam[k_hi,zr_i2:zr_i3])/\
                                     (V[k_low])*sgn*(dsm-i)/(dsm)
                for i in range(1,dsm):
                    F[k1,k2,zr_i3+i-1]=-c_*W[k1,k2]*dz_*ur2*\
                                     np.mean(lam[k_hi,zr_i2:zr_i3])/\
                                     (V[k_low])*sgn*(dsm-i)/(dsm)
                continue
            # level by level case
            # FIXME we have an issue here; mass balance 
            # (ie ignoring T and rtm just looking at mass of air)
            # is off on a layer by layer basis
            test=0
            for z in range(nz_):
                if z<=(zr_i-dsm):
                    ur = np.abs(T[k1,z]-T[k2,z])/T0*9.81**(.5)*l**(.5)
                    F[k1,k2,z]=c_*W[k1,k2]*dz_*ur*lam[k_low,z]/V[k_hi]*sgn
                    F_sum[k1,k2] = F_sum[k1,k2]+F[k1,k2,z]
                elif z<zr_i:
                    ur = np.abs(T[k1,z]-T[k2,z])/T0*9.81**(.5)*l**(.5)
                    F[k1,k2,z]=c_*W[k1,k2]*dz_*ur*lam[k_low,z]*(zr_i-z)/dsm/V[k_hi]*sgn
                elif z==zr_i:
                    continue
                elif z<(zr_i+dsm):
                    F[k1,k2,z]=-F[k1,k2,zr_i-(z-zr_i)]
                else:
                    test=test+1
                    F[k1,k2,z]=-1/denom*F_sum[k1,k2]
                    F_sumout[k1,k2] = F_sumout[k1,k2]+F[k1,k2,z]
    print(F[0,1,0])
    print(F[0,1,-1])
    print()
    return F
                    
    #FIXME consider virtual potential temperature vs temperature!!
'''

#### COMPUTE CIRCULATION FLUX ####
def circ_flux(W,T,lam,c_,H,V,vpt,k=k,T0=T0,l=l_het,dzi=dzi,dz_=dz):
    """
    calculate the circulating flux and assign its height/thickness
    
    Parameters
    ----------
    W  : int(k,k)
        Width of connection between two columns
    T  : float(k,nz)
        The potential temperature
    lam: float(k,nz)
        The species being fluxed; usually temperature or moisture
    c_ : float
        Constant flux parameter
    H  : float(k)
        Area average surface sensible heat flux
    V  : float(k)
        Volume of one layer
    vpt: float(k,nz)
        Virtual potential temperature
    k  : int
        Number of columns
    T0 : float
        Reference Temperature
    l  : float
        Lengthscale of heterogeneity
    dzi: int
        Thickness of circulation in gridspace
    dz_: int
        Thickness of one layer in meters

    Returns
    -------
    float(k,k,nz)
        Flux from k1 to k2. (+) is net flux out of k1 
    
    """
    
    # add the filter to make it smoother (beta)
    a = 1.5
    b = 1.5
    adj = beta.pdf(np.linspace(beta.ppf(.001,a,b),beta.ppf(.999,a,b),dzi),a,b)
    adj =adj/np.sum(adj)*dzi

    F = np.zeros((k,k,vpt.shape[1]))
    for k1 in range(k):
        for k2 in range(k):
            if k1==k2:
                continue
            sgn = (H[k2]-H[k1])/np.abs(H[k2]-H[k1])
            if sgn<0:
                k_low=k2
                k_hi =k1
            else:
                k_low=k1
                k_hi =k2

            # identify circulation height #FIXME only works (well) for k=2
            minz = int(np.where(vpt[k_hi,:]-vpt[k_low,:]<0)[0][0])
            
            T1m = np.mean(T[k1,0:dzi])
            T2m = np.mean(T[k2,0:dzi])
            ur = np.abs(T1m-T2m)/T0*9.81**(.5)*l**(.5)
            
            F[k1,k2,0:dzi] = c_*W[k1,k2]*dz_*ur*np.mean(lam[k_low,0:dzi])\
                             /(V[k_hi])*sgn*adj
            F[k1,k2,minz:minz+dzi] = -c_*W[k1,k2]*dz_*ur*\
                           np.mean(lam[k_hi,minz:minz+dzi])/(V[k_hi])*sgn*adj+\
                           F[k1,k2,minz:minz+dzi]
    print(F)
    print(V)
    print(W)
    print()
    return F


#### FIND NEAREST VALUE IN AN ARRAY ####
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
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
                data[x]=np.zeros((t,nz_))
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
    fp.write('! The vertical coordinate is entered in the first column as' +\
             ' height (z[m]) or pressure (Press[Pa]).\n')
    fp.write('! Temperature can be entered as liquid water potential '+\
             'temperature (thlm_f[K\s]), potential temperature (thm_f[K\s])'+\
             ', or absolute temperature(T_f[K\s]).\n')
    fp.write('! Prescribed, time-independent vertical velocity can be '+\
             'entered as velocity (w[m\s]) or pressure velocity '+\
             '(omega[Pa\s] or omega[mb\hr])\n')
    fp.write('Press[Pa]       T_f[K\s]        rtm_f[kg\kg\s]   '+\
            'um_ref[m\s]   vm_ref[m\s]      um_f[m\s^2]     '+\
            'vm_f[m\s^2]     omega[Pa\s]          ug[m\s]         vg[m\s]\n')
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


##############################################################################

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
    elif line[0:6]=='deltaz':
        dzline = 'deltaz  = '+str(dz)+'   ! Distance between grid levels on'+\
                 ' evenly-spaced grid.      [m]'
        lines.append(dzline+'\n')
    elif line[0:6]=='zm_top':
        zmline = 'zm_top = '+str(zmax)+' ! Maximum Altitude of '+\
                 'highest momentum level on any grid. [m]'
        lines.append(zmline+'\n')
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

# Fix tunable parameters, if necessary
if c_7 == -1:
    pass
else:
    lines=[]
    fp = open(blank_run+'input/configurable_model_flags.in','r')
    for line in fp:
        if line[0:8] == 'l_use_C7':
            lines.append('l_use_C7_Richardson          = .false.,\n')
        else:
            lines.append(line)
    fp.close()
    fp = open(m_dir+'/c_1/input/configurable_model_flags.in','w')
    for line in lines:
        fp.write(line)
    fp.close()
    
    lines=[]
    fp = open(blank_run+'input/tunable_parameters.in','r')
    for line in fp:
        if line[0:3] == 'C7 ':
            lines.append('C7          = '+str(c_7)+'00 ! Low Skewness in C7 '+\
                         'Skewness Function. Units [-]\n')
        elif line[0:3] == 'C7b':
            lines.append('C7b         = '+str(c_7)+'00 ! High Skewness in C7 '+\
                         'Skewness Function. Units [-]\n')
        else:
            lines.append(line)
    fp.close()
    fp = open(m_dir+'/c_1/input/tunable_parameters.in','w')
    for line in lines:
        fp.write(line)
    fp.close()


for i in list(range(2,k+1)):
    shutil.copytree(m_dir+'/c_1',m_dir+'/c_'+str(i))
    bindir = w_dir+'/k_'+str(k)+'/c_'+str(i)+'/bin/'
    os.remove(bindir+'clubb_standalone')
    os.symlink(cbin_dir+'clubb_standalone',bindir+'clubb_standalone')


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


print('WRITE SURFACE FILES',flush=True)

# write initial surface files using Nate's script
for i in list(range(k)):
    dir_old = os.getcwd()
    os.chdir(w_dir+'/k_'+str(k)+'/c_'+str(i+1)+'/other_scripts')
    subprocess.run('python create_arm_data_cpl.py',shell=True)
    os.chdir(dir_old)


# extend arm_forcing.in
nz_forcing = 37 #number of levels in original forcing file
for j in list(range(1,k+1)):
    frc_path = w_dir+'/k_'+str(k)+'/c_'+str(j)+\
                   '/input/case_setups/arm_forcings.in'
    data = read_forcings(frc_path)
    tmin = data['time'][0]
    tmax = data['time'][-1]
    nt_frc = (tmax-tmin)/(n_rest*delta_t)+1
    times_frc = np.linspace(int(tmin),int(tmax),int(nt_frc))
    data2 = {}
    for d in data.keys():
        if d == 'time':
            data2[d]=times_frc
        else:
            data2[d]=np.zeros((int(nt_frc),nz_forcing))
            for l in range(nz_forcing):
                data2[d][:,l]=np.interp(times_frc,data['time'],data[d][:,l])
    write_forcings(data2,frc_path)


# ------- #
# CLUSTER #
# ------- #
# cluster based on sensible heat and then adjust the other surface var
print('CLUSTERING',flush=True)

k_masks=np.zeros((nt,nx*nx))
for i in list(range(k)):
    Hv2 = Hv[np.mean(Hv,axis=1)>20,:]
    Hvn   = normalize(np.mean(Hv2[:,:],axis=0))
    X = np.reshape(Hvn,(-1,1))
    model = sklearn.cluster.KMeans(n_clusters=k)
    model.fit(X)
    Y = model.predict(X)
    for t in range(nt):
        k_masks[t,:] = Y[:]



# ------------------- #
# WRITE SURFACE FILES #
# ------------------- #
# transfer clustered array to surface files
print('OVERWRITE SURFACE FILES',flush=True)

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


# --------------- #
# OUTPUT CLUSTERS #
# --------------- #
print('OUTPUT CLUSTERS and COMPUTE CONNECTIONS')
fp=nc.Dataset(w_dir+'/k_'+str(k)+'/clusters.nc','w')
fp.createDimension('t',size=nt)
fp.createDimension('lat',size=nx)
fp.createDimension('lon',size=nx)
fp.createDimension('clst',size=k)
fp.createVariable('W','d',dimensions=('clst','clst'))
fp.createVariable('cluster','d',dimensions=('t','lon','lat'))
fp.createVariable('tskin','d',dimensions=('t','lon','lat'))
fp.createVariable('H','d',dimensions=('t','lon','lat'))
fp.createVariable('LE','d',dimensions=('t','lon','lat'))
fp.createVariable('frac','d',dimensions=('clst'))
fp.createVariable('H_clst','d',dimensions=('t','clst'))

clst_2d = np.reshape(k_masks[:,:],(nt,nx,nx))
fp['tskin'][:]=Tsfcg[:]
fp['cluster'][:]=clst_2d[:]
fp['H'][:]=Hg[:]
fp['LE'][:]=Lg[:]

# compute fraction that each cluster occupies
clst_frac = np.zeros((k,))
for i in range(k):
    clst_frac[i]=np.sum(k_masks[:]==i)/k_masks.size
fp['frac'][:]=np.array(clst_frac)[:]

# compute the width of the connecting surface 
W = np.zeros((k,k))
for i in range(nx):
    for j in range(nx):
        cl = int(clst_2d[0,i,j])
        if i>0:
            cl_2 = int(clst_2d[0,i-1,j])
            W[cl,cl_2]=W[cl,cl_2]+dx
        if i<(nx-1):
            cl_2 = int(clst_2d[0,i+1,j])
            W[cl,cl_2]=W[cl,cl_2]+dx
        if j>0:
            cl_2 = int(clst_2d[0,i,j-1])
            W[cl,cl_2]=W[cl,cl_2]+dx
        if j<(nx-1):
            cl_2 = int(clst_2d[0,i,j+1])
            W[cl,cl_2]=W[cl,cl_2]+dx
fp['W'][:] = W[:]

# compute the sensible heat flux of each cluster
H_clst = np.zeros((nt,k))
for i in range(k):
    for t in range(nt):
        H_clst[t,i]=np.mean(fp['H'][t,:,:][fp['cluster'][t,:,:]==i])
fp['H_clst'][:]=H_clst[:]

print(np.mean(H_clst,axis=0))

# compute the volume of each cluster
V = clst_frac*dz*nx*nx*dx*dx

#sys.exit()

# --------- #
# CORE LOOP #
# --------- #
tlist = list(range(int(t_init),int(t_final),int(round(delta_t*n_rest))))
H2 = np.zeros((len(tlist),k))
tss = []
for t in range(nt):
    tss.append(t_init+t*dt)
for i in range(k):
    H2[:,i]=np.mean(H_clst[:,i])
    # uncomment below for variable hot and warm patches
    #H2[:,i]=np.interp(tlist,tss,H_clst[:,i])

# FIRST TIMESTEP 
for j in list(range(1,k+1)):
    # run the script
    cbasedir= w_dir+'/k_'+str(k)+'/c_'+str(j)
    rfile = cbasedir+'/run_scripts/run_scm.bash'
    print(rfile)
    cmd = rfile+' arm -p '+cbasedir+'/input/tunable_parameters.in'+' >'+'log'+str(j)
    subprocess.run(cmd,shell=True,stdin=subprocess.DEVNULL)
    print()
    print(str(j) + ' IS COMPLETE ',flush=True)
    os.mkdir(w_dir+'/k_'+str(k)+'/c_'+str(j)+'/restart')
    os.mkdir(w_dir+'/k_'+str(k)+'/c_'+str(j)+'/output/old')

frcs = read_forcings(m_dir+'/c_'+str(1)+'/input/case_setups/arm_forcings.in',nz_forcing)
p_in = frcs['Press[Pa]'][0,:]
t_in = frcs['time'][:]
frcs = {}

# OUTPUT RUN INFORMATION 
fp = open(w_dir+'/tw_run_param.txt','w')
fp.write('n_rest  = '+str(n_rest)+' # timesteps to next restart\n')
fp.write('delta_t = '+str(delta_t)+' # clubb timestep in seconds\n')
fp.write('dsmooth = '+str(dsmooth)+' # number of vertical levels to smooth forcing at switch\n')
fp.write('T0      = '+str(T0)+' # clubb reference temperature\n')
fp.write('c_t     = '+str(c_t)+' # flux factor for temperature\n')
fp.write('c_r     = '+str(c_r)+' # flux factor for moisture\n')
fp.write('l_het   = '+str(l_het)+' # lengthscale of heterogeneity\n')
fp.write('sfc_dir = '+sfc_dir+' # sfc directory\n')
fp.write('k       = '+str(k)+' # number of columns\n')
fp.write('nx      = '+str(nx)+' # number of gridcells in each direction\n')
fp.write('dx      = '+str(dx)+' # surface grid resolution\n')
fp.write('dz      = '+str(dz)+' # vertical grid resolution\n')
fp.write('zmax    = '+str(zmax)+' # maximum height in clubb\n')
fp.write('z_r     = '+str(z_r)+' # height where we switch sign of circ flux\n')
fp.write('dzc     = '+str(dzc)+' # thickness of circulation if no z_r used\n')
fp.write('stdate  = '+str(stdate)+' # start date\n')
fp.write('enddate = '+str(enddate)+' # end date\n')
fp.write('\nRan at: '+str(datetime.datetime.today()))
fp.close()

# EXTEND FORCING TO COVER FULL VERTICAL DOMAIN
for i in range(1,k+1):
    frcs_i = read_forcings(m_dir+'/c_'+str(i)+'/input/case_setups/arm_forcings.in',nz_forcing)
    zs = np.linspace(0,zmax,nz)
    pres_frc = frcs_i['Press[Pa]'][:]
    fp_output = nc.Dataset(m_dir+'/c_'+str(i)+'/output/arm_zt.nc','r')
    pres_newfrc = fp_output['p_in_Pa'][0,:,0,0]
    fp_output.close()
    data2 = {}
    for var in frcs_i.keys():
        if var == 'time':
            data2[var]=np.zeros((len(tlist)+1))
            data2[var][:]=frcs_i[var][:]
            continue
        data2[var]=np.zeros((len(tlist)+1,nz))
        if var == 'Press[Pa]':
            for t in range(len(tlist)):
                data2[var][t,:]=pres_newfrc[:]
            continue
        for t in range(len(tlist)):
            data2[var][t,:]=np.interp(pres_newfrc,pres_frc[2,:][::-1],frcs_i[var][t,:][::-1])
    write_forcings(data2,m_dir+'/c_'+str(i)+'/input/case_setups/arm_forcings.in')

# read in forcings
for i in range(1,k+1):
    frc_file = m_dir+'/c_'+str(i)+'/input/case_setups/arm_forcings.in'
    frcs[i] = read_forcings(frc_file,nz)
for i in range(1,len(tlist)): 
    t0 = tlist[i]
    t1 = t0+n_rest*delta_t
    tmps = np.zeros((k,nz))
    rtms = np.zeros((k,nz))
    thlm = np.zeros((k,nz))
    thvm = np.zeros((k,nz))

    # Load in the temperature and mixing ratio
    for j in range(1,k+1):
        fp = nc.Dataset(m_dir+'/c_'+str(j)+'/output/arm_zt.nc','r')
        tmps[j-1,:]=fp['T_in_K'][int(round(n_rest-1)),:,0,0]
        rtms[j-1,:]=fp['rtm'][int(round(n_rest-1)),:,0,0]
        thlm[j-1,:]=fp['thlm'][int(round(n_rest-1)),:,0,0]
        thvm[j-1,:]=fp['thvm'][int(round(n_rest-1)),:,0,0]
        fp.close()

    # check if there is surface heating; if no surface heating no flux
    doflux = True
    if frcs[1]['T_f[K\s]'][i,0]<0:
        doflux=False

    # compute and define fluxes
    if doflux:
        # OLD WAY
        #zri = int(np.floor((z_r)/dz)+1)
        #zri2 = int(np.floor((z_r2)/dz)+1)
        #zri3 = int(np.floor((z_r3)/dz)+1)
        #F_T = circ_flux(W,thlm,tmps,c_t,H2[i,:],V,zri,zri2,zri3)
        #F_r = circ_flux(W,thlm,rtms,c_r,H2[i,:],V,zri,zri2,zri3)
        
        F_T = circ_flux(W,thlm,tmps,c_t,H2[i,:],V,thvm)
        F_r = circ_flux(W,thlm,rtms,c_r,H2[i,:],V,thvm)

    for j in list(range(1,k+1)):
        frc_file = m_dir+'/c_'+str(j)+'/input/case_setups/arm_forcings.in'
        
        if doflux:
            # Change Flux to forcings
            frcs[j]['T_f[K\s]'][i,:]=frcs[j]['T_f[K\s]'][i,:]-\
                                   np.sum(F_T[j-1,:,:],axis=0)
            frcs[j]['rtm_f[kg\kg\s]'][i,:]=frcs[j]['rtm_f[kg\kg\s]'][i,:]-\
                                         np.sum(F_r[j-1,:,:],axis=0)
            write_forcings(frcs[j],frc_file)

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
                lines.append('time_restart = '+str(t0)+' \n')
            elif line[0:10]=='time_final':
                fln = 'time_final = '+str(t1).ljust(11)+\
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
        cbasedir= w_dir+'/k_'+str(k)+'/c_'+str(j)
        rfile = cbasedir+'/run_scripts/run_scm.bash'
        print(rfile)
        cmd = rfile+' arm -p '+cbasedir+'/input/tunable_parameters.in'
        
        print('STARTING RUN: '+str(j)+' at time '+str(t0)+' with process '+str(0),flush=True)
        subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
        print('RUN COMPLETE: '+str(j)+' at time '+str(t0)+' with process '+str(0),flush=True)

# ---------- #
# AGG OUTPUT #
# ---------- #
# temporally aggregate small output files
# then do traditional output aggregation


        

    
    

    
