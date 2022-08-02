#!/usr/bin/python
# Runs a crudely coupled model; 

# IMPORT 
import sys
import argparse
import subprocess
import numpy as np
import random
import os
import datetime
import sklearn.cluster
import netCDF4 as nc
import shutil
from scipy.stats import beta
from inspect import currentframe, getframeinfo

# ------------------------- #
# USER INPUTS and CONSTANTS #
# ------------------------- #

#### CONSTANTS ####
# run options
s_rest      = 300 # number of seconds between restarts
delta_t     = 30
dsmooth     = 5 # number of levels over which to switch sign of circ flux
T0          = 300 # reference temperature; matched to CLUBB
c_rat       = 1 # ratio of moisture to temperature factor
c_t         = 1 #.00001 # factor for temperature fluxes
c_ur        = .0001
hggt        = 2 # hours after start to use surface as basis for l_het
no_wind_frc = False # if true, will ignore forcings for wind velocity

# directories
sfc_dir   = '/stor/soteria/hydro/shared/sgp_surfaces/dx0100nx1000/'
clubb_dir = '/home/tsw35/tyche/clubb/' # directory of clubb run folders
blank_run = '/home/tsw35/soteria/clubb/clubb_scripts/run_/' # a clean run folder to copy
cbin_dir  = '/home/tsw35/soteria/software/CLUBB/bin/' # bin for clubb
met_dir   = '/home/tsw35/soteria/clubb/data/sgp60varanarap_2012-2019.nc' #forcing
tune_o    = 'X' #'/home/tsw35/soteria/clubb/clubb_scripts/tunable_param/g0.4c110.4c80.5' 
            # 'X' overwrite tunable parameters


#### DEFAULTS ####
k       = 2 # number of clusters
nx      = 1000 # number of gridcells in each direction
dx      = 100  # resolution of surface grid in meters
dz      = 40 # dz
zmax    = 12000
dzc     = 1000
nz      = int(np.floor(zmax/dz)+1)
stdate  = '2017-07-16T12:00:00.000'# start date iso format
enddate = '2017-07-17T03:00:00.000'# end date iso format
dirname = 'test_cpl6_nocpl'
l_het   = -2 # Length of heterog. in meters; -1: calculate, -2: || to wind
c_7     = -1 # default is non constant (-1), bouancy term, .3 to .8 is listed range
flux_on = 0 # turn on cirrculation flux (1) turn off (0)
wind_cr = 1 # turn on directional wind corrections to flux (1) turn off (0)
cut_off= .8 # cuttoff for correlation function heterogeneity
              # lowering increases the velocity of circulation

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
prs.add_argument('-r', action='store', dest='c_rat',type=float,default=c_rat)
prs.add_argument('-t', action='store', dest='c_t',type=float,default=c_t)
prs.add_argument('-v', action='store', dest='c_7',type=float,default=c_7)
prs.add_argument('-f', action='store', dest='flux',type=int,default=flux_on)
prs.add_argument('-w', action='store', dest='wc',type=int,default=wind_cr)
prs.add_argument('-y', action='store', dest='dt',type=float,default=delta_t)
prs.add_argument('-o', action='store', dest='tun',default=tune_o)
prs.add_argument('-u', action='store', dest='cut',type=float,default=cut_off)
prs.add_argument('-a', action='store', dest='atm',default='')

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
c_rat   = args.c_rat
c_t     = args.c_t
c_r     = c_t*c_rat
c_7     = args.c_7
flux_on = args.flux==1
wind_cr = args.wc==1
delta_t = args.dt
tune_o  = args.tun #overwrite tunable param w/ other file
cut_off = args.cut
atm_dir = args.atm

#### EXTRAS ####
dzi = int(np.floor(dzc/dz))
n_rest = int(np.floor(s_rest/delta_t)) # number of timesteps between restart
if k == 1:
    flux_on = False


##############################################################################

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


#### COMPUTE CIRCULATION FLUX ####
def circ_flux(W,T,lam,c_,H,V,vpt,k=k,T0=T0,l=l_het,dzi=dzi,dz_=dz,w_=[0,0],cd_=np.zeros((k,k,2))):
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
    w_ : float(2)
        Mean wind velocity of lowest 1000 m of atmosphere
    cd_: float(k,k,2)
        Directionality of connection between two columns; x,y

    Returns
    -------
    float(k,k,nz)
        Flux from k1 to k2. (+) is net flux out of k1 
    float(k,nz)
        Vertical velocity between each layer; (+) is up; m/s
    float
        The velocity of the circulation
    float
        The base height of the circulation
    """
    ## START WITH HORIZONTAL FLUX ## 
    # add the filter to make it smoother (beta)
    a = 1.5
    b = 1.5
    print('c_:%f'%c_)
    adj = beta.pdf(np.linspace(beta.ppf(.001,a,b),beta.ppf(.999,a,b),dzi),a,b)
    adj =adj/np.sum(adj)*dzi
    
    F_v = np.zeros((k,k,vpt.shape[1])) #volumetric flux of air
    F = np.zeros((k,k,vpt.shape[1]))
    F_z = np.zeros((k,vpt.shape[1]))

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
            try:
                minz = int(np.where(vpt[k_hi,:]-vpt[k_low,:]<0)[0][0])
            except:
                return F,F_z,0,0,float('nan'),float('nan')
            print('minz: %f'%minz)
            print('k_hi: %f'%k_hi)
            print('k_low: %f'%k_low)
            T1m = np.mean(T[k1,0:dzi])
            T2m = np.mean(T[k2,0:dzi])
            
            ur = c_ur*np.abs(T1m-T2m)/T0*9.81**(.5)*l**(.5)
            
            # adjust for velocity of mean wind
            if (w_[0]==0) and (w_[1]==0):
                pass
            else:
                # have mean wind reduce circulation velocity
                '''
                circ_dir=cd_[k_low,k_hi,:]
                flow=circ_dir/np.sqrt(circ_dir[0]**2+circ_dir[1]**2)*ur
                adj_flow = np.abs(flow)-np.abs(w_)*.66
                if adj_flow[0]<0:
                    adj_flow[0]=0
                if adj_flow[1]<0:
                    adj_flow[1]=0
                ur=np.sqrt(adj_flow[0]**2+adj_flow[1]**2)
                '''

                # take circulation velocity component perp to wind
                circ_dir=cd_[k_low,k_hi,:]
                flow=circ_dir/np.sqrt(circ_dir[0]**2+circ_dir[1]**2)*ur
                adj_flow = flow-np.dot(np.dot(flow,w_)/np.dot(w_,w_),w_)
                ur=np.sqrt(adj_flow[0]**2+adj_flow[1]**2)
            
            
            # compute the flux
            F_v[k1,k2,0:dzi]=W[k1,k2]*dz_*ur*sgn*adj
            F_v[k1,k2,minz:minz+dzi]=-W[k1,k2]*dz_*ur*sgn*adj

            F[k1,k2,0:dzi] = c_*np.mean(lam[k_low,0:dzi])*\
                             F_v[k1,k2,0:dzi]/(V[k1])
            F[k1,k2,minz:minz+dzi] = c_*np.mean(lam[k_hi,minz:minz+dzi])/\
                    (V[k1])*F_v[k1,k2,minz:minz+dzi]+F[k1,k2,minz:minz+dzi]
    xfac=np.abs(np.mean(F[0,1,0:dzi]))
    
    ## VERTICAL FLUX ##
    # determined by pure volumetric flux balance
    A = V/dz_ # area of a layer
    
    
    for k_ev in range(k):
        out_upv=0 # volume flux leaving to (+) or entering from (-) layer abve
        for z in range(minz+dzi):
            in_downv= out_upv # flx entering frm (+) leaving to (-) layr below
            in_hv=0 # volumetric only
            for k2 in range(k):
                if k_ev==k2:
                    continue
                in_hv=in_hv+F_v[k2,k_ev,z]
            
            # balance the flux from layer below and other column
            
            out_upv = (in_hv+in_downv)
            
            F_z[k_ev,z]=out_upv/A[k_ev]
    
    urz = np.max(np.abs(F_z),axis=1)
    print()
    print('ur: %f'%ur)
    print('urz: %f,%f'%(urz[0],urz[1]))
    print()
    return F,F_z,ur,urz,minz,xfac


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


#### CONVERT WS [M/S] TO OMEGA [Pa/S] ####
def w_to_omega(ws,press,tmp):
    '''
        Parameters

        ws     : N,1 array, vertical windspeed in m/s
        press  : N,1 array, pressure in Pa
        tmp    : N,1 array, temperature
        
    '''
    rgas=287.058
    g=9.8067
    rho = press/(rgas*tmp)
    omega=-ws*rho*g
    return omega



#### CALCULATE LENGTHSCALE OF HETEROGENEITY ####
def estimate_l_het(l_het,Hg_,cut=.25,samp=10000):
    '''
        Parameters

        l_het  : -1: compute full lengthscale, -2: compute parallel to wind
        cut    : percentage cuttoff for lengthscale
        Hg     : a N by N grid of the sensible heat fluxes
        samp   : how many points to sample when constructing cov matrix
                 sample size is 10 times this when parallel to wind
    '''

    # CASE 0: Regular constant lengthscale
    if l_het>-1:
        return l_het

    # Common Setup for Cases 1 and 2
    a_,b_ = Hg_.shape
    r_H = np.zeros((a_,b_))
    c_H = np.zeros((a_,b_))
    for i in range(a_):
        for j in range(b_):
            r_H[i,j]=i*dx
            c_H[i,j]=j*dx
    H_flat = Hg_.flatten()
    r_flat = r_H.flatten()
    c_flat = c_H.flatten()


    # CASE 1: Full heterogneeity 
    if l_het == -1:

        idx = np.random.choice(len(H_flat),size=samp,replace=False)
        H_sb = H_flat[idx]
        r_Hsb = r_flat[idx]
        c_Hsb = c_flat[idx]
        mu=np.mean(H_sb)
        a = (H_sb[:,np.newaxis].T - mu)*(H_sb[:,np.newaxis]-mu)
        h = ((r_Hsb[:,np.newaxis] - r_Hsb.T)**2 + \
            (c_Hsb[:,np.newaxis] - c_Hsb.T)**2)**0.5
        Qf = a.flatten()
        hf = h.flatten()

        bins = np.linspace(0,50000,51)
        means=np.zeros((len(bins)-1,))
        for i in range(len(bins)-1):
            means[i]=np.mean(Qf[(hf>bins[i])&(hf<bins[i+1])])

        l_het_ = bins[0:-1][means<=(.25*means[0])][0]
        
        print('Non Dimensional Lengthscale of Heterogeneity')
        print(str(l_het_))
        print()
        

    # CASE 2: In Direction of Mean Wind
    elif l_het == -2:
        
        # import mean wind information
        sound_path = w_dir+'/k_'+str(k)+'/c_'+str(k)+\
                   '/input/case_setups/arm_sounding.in'
        fp=open(sound_path,'r')
        u_=0
        v_=0
        for line in fp:
            if line[0]=='!':
                continue
            if line[0]=='P':
                continue
            if line[0]=='z':
                continue
            linesp = line.split(' ')
            print(linesp)
            ll=[]
            for l_ in linesp:
                if not l_=='':
                    ll.append(l_)
            u_=float(ll[3])
            v_=float(ll[4])
            break
        fp.close()
        
        # add random mean wind to avoid code issues under synthetic cases
        u_=u_+random.uniform(-.0001*u_,.0001*u_)
        v_=v_+random.uniform(-.0001*v_,.0001*v_)
        
        # normalize mean wind
        u_p = u_/(u_**2+v_**2)**(1/2)
        v_p = v_/(u_**2+v_**2)**(1/2)
        
        # select half the points; compute mu
        idx = np.random.choice(len(H_flat),size=int(round(len(H_flat)/2)),replace=False)
        mu = np.mean(H_flat[idx])

        # x heterogeneity
        if np.abs(u_p)>(1/np.sqrt(2)):
            Qf=[]
            hf=[]
            for i in range(a_):
                idx_i = idx[r_flat[idx]==i]
                H_sbi = H_flat[idx_i]
                c_Hsbi= c_flat[idx_i]
                a = (H_sbi[:,np.newaxis].T - mu)*(H_sbi[:,np.newaxis]-mu)
                h = ((c_Hsbi[:,np.newaxis] - c_Hsbi.T)**2)**0.5
                Qf.extend(a.flatten())
                hf.extend(h.flatten())
            bins = np.linspace(0,100000,101)
            means=np.zeros((len(bins)-1,))
            Qf=np.array(Qf)
            hf=np.array(hf)
            for j in range(len(bins)-1):
                means[j]=np.mean(Qf[(hf>bins[j])&(hf<bins[j+1])])
            try:
                l_het_a = bins[0:-1][means<=(.25*means[0])][0]
            except:
                l_het_a = 100000
            a_i = 1
            a_j = 0
                
        
        # xy heterogeneity
        if (u_p*v_p)>0:
            Qf=[]
            hf=[]
            for i in range(-round(min(a_,b_)/2),round(min(a_,b_)/2)):
                idx_i = idx[(c_flat[idx]+i)==r_flat[idx]]
                H_sbi = H_flat[idx_i]
                c_Hsbi= c_flat[idx_i]
                a = (H_sbi[:,np.newaxis].T - mu)*(H_sbi[:,np.newaxis]-mu)
                h = ((c_Hsbi[:,np.newaxis] - c_Hsbi.T)**2)**0.5
                Qf.extend(a.flatten())
                hf.extend(h.flatten())
            Qf=np.array(Qf)
            hf=np.array(hf)
            bins = np.linspace(0,100000,101)
            means=np.zeros((len(bins)-1,))
            for j in range(len(bins)-1):
                means[j]=np.mean(Qf[(hf>bins[j])&(hf<bins[j+1])])
            try:
                l_het_b = bins[0:-1][means<=(.25*means[0])][0]
            except:
                l_het_b = 50000*np.sqrt(2)
                
            # i and j components of xy vector
            b_i=1/np.sqrt(2)
            b_j=1/np.sqrt(2)
        
        # x-y heterogeneity
        if (u_p*v_p)<0:
            Qf=[]
            hf=[]
            for i in range(round(min(a_,b_)/2),round(min(a_,b_)/2)*3):
                idx_i = idx[(c_flat[idx]+i)==r_flat[idx]]
                H_sbi = H_flat[idx_i]
                c_Hsbi= c_flat[idx_i]
                a = (H_sbi[:,np.newaxis].T - mu)*(H_sbi[:,np.newaxis]-mu)
                h = ((c_Hsbi[:,np.newaxis] - c_Hsbi.T)**2)**0.5
                Qf.extend(a.flatten())
                hf.extend(h.flatten())
            Qf=np.array(Qf)
            hf=np.array(hf)
            bins = np.linspace(0,100000,101)
            means=np.zeros((len(bins)-1,))
            for j in range(len(bins)-1):
                means[j]=np.mean(Qf[(hf>bins[j])&(hf<bins[j+1])])
            try:
                l_het_b = bins[0:-1][means<=(.25*means[0])][0]
            except:
                l_het_b = 50000*np.sqrt(2)
                
            # i and j components of -xy vector
            b_i=1/np.sqrt(2)
            b_j=-1/np.sqrt(2)
        

        # y heterogeneity 
        if np.abs(v_p)>(1/np.sqrt(2)):
            Qf=[]
            hf=[]
            for i in range(b_):
                idx_i = idx[c_flat[idx]==i]
                H_sbi = H_flat[idx_i]
                r_Hsbi= r_flat[idx_i]
                a = (H_sbi[:,np.newaxis].T - mu)*(H_sbi[:,np.newaxis]-mu)
                h = ((r_Hsbi[:,np.newaxis] - r_Hsbi.T)**2)**(.5)
                Qf.extend(a.flatten())
                hf.extend(h.flatten())

            bins = np.linspace(0,100000,101)
            means=np.zeros((len(bins)-1,))
            Qf=np.array(Qf)
            hf=np.array(hf)
            for j in range(len(bins)-1):
                means[j]=np.mean(Qf[(hf>bins[j])&(hf<bins[j+1])])
            try:
                l_het_a = bins[0:-1][means<=(.25*means[0])][0]
            except:
                l_het_a = 100000
            
            # i and j components of y vector
            a_i=0
            a_j=1
                               
        # final computation; two cases depending on a_i
        if a_i==1:
            beta  = (v_p*a_i-u_p*a_j)/(b_j-b_i*a_j)
            alpha = (u_p-beta*b_i)/a_i
        else:
            beta  = (u_p*a_j-v_p*a_i)/(b_i-b_j*a_i)
            alpha = (v_p-beta*b_j)/a_j 
        
        l_het_=alpha*l_het_a+beta*l_het_b
        if l_het_<0:
            print('ERROR POSSIBLE IN HETEROGENEITY')
            print('ERROR POSSIBLE IN HETEROGENEITY')
            print('TAG:NEGATIVE_HET')
            l_het_=np.abs(l_het_)
        
        print('Heterogeneity in direction x:'+str(u_p)+'  y:'+str(v_p))
        print(l_het_)
        print('')
    return l_het_



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
    elif line[0:7]=='dt_main':
        dtln= 'dt_main = '+str(delta_t)+' ! Model timestep [s]'
        lines.append(dtln+'\n')
    elif line[0:6]=='dt_rad':
        dtln= 'dt_rad = '+str(delta_t)+' ! Model timestep [s]'
        lines.append(dtln+'\n')
    elif line[0:10]=='stats_tout':
        dtln= 'stats_tout = '+str(max(60,delta_t))+' ! Statistical Sampling Timestep [s]'
        lines.append(dtln+'\n')
    elif line[0:11]=='stats_tsamp':
        dtln= 'stats_tsamp = '+str(delta_t)+' ! Statistical Sampling Timestep [s]'
        lines.append(dtln+'\n')
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
if not (tune_o == 'X'):
    subprocess.run('cp '+tune_o+' '+m_dir+'/c_1/input/tunable_parameters.in',shell=True)
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


# Read in the surface data
# assume data is hourly 
print('... read in surface data',flush=True)
nt     = int((endt-stdt).seconds/3600)+1
Hg,Hv   = read_sfc_data('sh',nt,stdt)
Lg,Lv   = read_sfc_data('lh',nt,stdt)
lwg,lwv  = read_sfc_data('lw',nt,stdt)
#lon_g,lonv = read_sfc_data('',1,stdt,sfc_dir+'jsslongrid_02')
#lat_g,latv = read_sfc_data('',1,stdt,sfc_dir+'jsslatgrid_02')
Tsfcv = (lwv/(5.67*10**(-8)))**(1/4)
Tsfcg = (lwg/(5.67*10**(-8)))**(1/4)

print('WRITE SURFACE FILES',flush=True)

# write initial surface files using Nate's script
for i in list(range(k)):
    dir_old = os.getcwd()
    os.chdir(w_dir+'/k_'+str(k)+'/c_'+str(i+1)+'/other_scripts')
    if atm_dir=='':
        subprocess.run('python create_arm_data_cpl.py',shell=True)
    else:
        subprocess.run('python create_arm_data_cpl.py -a '+atm_dir,shell=True)
    os.chdir(dir_old)

# pull in the sw down
t0_var = datetime.datetime(2012,5,1,0,0)
fp_sw=nc.Dataset(met_dir,'r')
swt_init =int(((stdt-t0_var).total_seconds()/3600+1)) #initial time to read
swt_final=int(((endt-t0_var).total_seconds()/3600+1)) #final time to read
sw_dwn =fp_sw['sw_dn_srf'][swt_init:swt_final]

print('...Calculate Lengthscale of Heterogeneity',flush=True)
# Calculate Lengthscale of Heterogeneity
Hgg = Hg[hggt,:,:]
l_het = estimate_l_het(l_het,Hgg,cut=cut_off)

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
        elif ((d == 'vm_ref[m\s]') or (d == 'omega[Pa\s]') or (d == 'um_ref[m\s]')) and no_wind_frc:
            data2[d]=np.ones((int(nt_frc),nz_forcing))*-999.9
        else:
            data2[d]=np.zeros((int(nt_frc),nz_forcing))
            for l in range(nz_forcing):
                data2[d][:,l]=np.interp(times_frc,data['time'],data[d][:,l])
    
    # save full structure aka "original" and other
    write_forcings(data2,frc_path)
    write_forcings(data2,w_dir+'/arm_forcings_o.in')


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
cdirect = np.zeros((k,k,2))
for i in range(nx):
    for j in range(nx):
        cl = int(clst_2d[0,i,j])
        if i>0:
            cl_2 = int(clst_2d[0,i-1,j])
            W[cl,cl_2]=W[cl,cl_2]+dx
            cdirect[cl,cl_2,1]=cdirect[cl,cl_2,1]+1
        if i<(nx-1):
            cl_2 = int(clst_2d[0,i+1,j])
            W[cl,cl_2]=W[cl,cl_2]+dx
            cdirect[cl,cl_2,1]=cdirect[cl,cl_2,1]-1
        if j>0:
            cl_2 = int(clst_2d[0,i,j-1])
            W[cl,cl_2]=W[cl,cl_2]+dx
            cdirect[cl,cl_2,0]=cdirect[cl,cl_2,0]-1
        if j<(nx-1):
            cl_2 = int(clst_2d[0,i,j+1])
            W[cl,cl_2]=W[cl,cl_2]+dx
            cdirect[cl,cl_2,0]=cdirect[cl,cl_2,0]+1
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

##############################################################################

# --------- #
# CORE LOOP #
# --------- #
# timestep list
tlist = list(range(int(t_init),int(t_final),int(round(delta_t*n_rest))))

#interpolate shortwave down
tlistsw=list(range(int(t_init),int(t_final),3600))
sw_dwn_interp = np.interp(tlist,tlistsw,sw_dwn)
print(sw_dwn_interp)

# create an array of flux activation w/ 1 active, 0 inactive, between is decay
min_rad = 250
fluxbool=np.zeros((len(sw_dwn_interp),))
fluxbool[sw_dwn_interp>min_rad]=1
if flux_on:
    # change boolean to decay setup
    doflux=[0]
    for i in range(1,len(fluxbool)):
        prev=doflux[i-1]
        if (prev<.99)&(fluxbool[i]==1):
            doflux.append(prev+.1667)
        elif (prev>.01)&(fluxbool[i]==0):
            doflux.append(prev-.1667)
        elif fluxbool[i]==0:
            doflux.append(0)
        elif fluxbool[i]==1:
            doflux.append(1)
        else:
            print('ERROR -- FLUXBOOL NOT 0 or 1')
else:
    doflux=np.zeros((len(sw_dwn_interp),))
# patches 
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

# Initialize circulation output variables
u_r = np.ones((len(tlist),))*-1
u_rz=np.ones((2,len(tlist)))*-1
z_circ = np.ones((len(tlist),))*-1
fpow = np.ones((len(tlist),))*-1

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
    um   = [0,0]

    try:
        # Load in the temperature and mixing ratio
        n_rest_out=int(np.floor(s_rest/max(60,delta_t)))
        for j in range(1,k+1):
            fp = nc.Dataset(m_dir+'/c_'+str(j)+'/output/arm_zt.nc','r')
            tmps[j-1,:]=fp['T_in_K'][int(round(n_rest_out-1)),:,0,0]
            rtms[j-1,:]=fp['rtm'][int(round(n_rest_out-1)),:,0,0]
            thlm[j-1,:]=fp['thlm'][int(round(n_rest_out-1)),:,0,0]
            thvm[j-1,:]=fp['thvm'][int(round(n_rest_out-1)),:,0,0]
        
            # load in windspeed
            if wind_cr:
                um[0]=np.mean(fp['um'][int(round(n_rest_out-1)),0:25,0,0])
                um[1]=np.mean(fp['vm'][int(round(n_rest_out-1)),0:25,0,0])

            fp.close()
    except KeyboardInterrupt:
        sys.exit()
        pass
    except Exception as e:
        print(e)
        frameinfo = getframeinfo(currentframe())
        print(frameinfo.filename, frameinfo.lineno)
        break
    
    #depricated flux controller (doflux) was once here

    # compute and define fluxes
    if doflux[i]>.01:
        F_T,F_wz,u_r[i],u_rz[:,i],z_circ[i],fpow[i] = circ_flux(W,thlm,tmps,c_t*doflux[i],H2[i,:],V,thvm,\
                                         l=l_het,w_=um,cd_=cdirect)
        F_r = circ_flux(W,thlm,rtms,c_r*doflux[i],H2[i,:],V,thvm,l=l_het,\
                        w_=um,cd_=cdirect)[0]

    for j in list(range(1,k+1)):
        frc_file = m_dir+'/c_'+str(j)+'/input/case_setups/arm_forcings.in'
       
        # Change fluxes to forcings
        if doflux[i]>.01:
            frcs[j]['T_f[K\s]'][i,:]=frcs[j]['T_f[K\s]'][i,:]-\
                                   np.sum(F_T[j-1,:,:],axis=0)
            frcs[j]['rtm_f[kg\kg\s]'][i,:]=frcs[j]['rtm_f[kg\kg\s]'][i,:]-\
                                         np.sum(F_r[j-1,:,:],axis=0)
            # change from m/s to Pa/s
            omega=w_to_omega(F_wz[j-1,:],frcs[j]['Press[Pa]'][i,:],tmps[j-1,:])
            frcs[j]['omega[Pa\s]'][i,:]=frcs[j]['omega[Pa\s]'][i,:]-omega
        
        # Pull out required forcing timesteps
        frcs_sm = {}
        for var in frcs[j].keys():
            if var == 'time':
                frcs_sm[var]=frcs[j][var][i:i+2]
            else:
                frcs_sm[var]=frcs[j][var][i:i+2,:]
        write_forcings(frcs_sm,frc_file)

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
        try:
            subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
        except Exception as e:
            print('ERROR RUN FAILED')
            print(e)
            break
        print('RUN COMPLETE: '+str(j)+' at time '+str(t0)+' with process '+str(0),flush=True)

#### OUTPUT ####
# Output Final forcings
for j in list(range(1,k+1)):
    write_forcings(frcs[j],w_dir+'/k_'+str(k)+'/c_'+str(j)+'/arm_forcings_f.in')

# output circulation characteristics
fp=nc.Dataset(w_dir+'/k_'+str(k)+'/clusters.nc','r+')
fp.createDimension('xy',size=2)
fp.createDimension('t_model',size=len(tlist))
fp.createVariable('t_model','d',dimensions=('t_model'))
fp.createVariable('u_r','d',dimensions=('t_model'))
fp.createVariable('z_circ','d',dimensions=('t_model'))
fp.createVariable('flux_power','d',dimensions=('t_model'))
fp.createVariable('flow_dir','d',dimensions=('clst','clst','xy'))
fp.createVariable('u_r_z','d',dimensions=('clst','t_model'))
fp['u_r'][:]=u_r[:]
fp['flux_power'][:]=fpow[:]
fp['z_circ'][:]=z_circ[:]*dz
fp['flow_dir'][:]=cdirect[:]
fp['u_r_z'][:]=u_rz[:]

        

    
    

    

