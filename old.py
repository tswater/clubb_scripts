#!/usr/bin/python
# Runs a crudely coupled model; 

# IMPORT 
import sys
import rasterio
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
decay_t     = 60*120 # number of seconds to full circulation power
delta_t     = 30
T0          = 300 # reference temperature; matched to CLUBB
hggt        = 6 # hours after start to use surface as basis for l_het
vert_circ   = True # if true, will include vertical circulation
no_tq_frc   = False # if true, will set thermal and moisture forcings to 0


# directories
sfc_dir   = '/home/tsw35/soteria/clubb/data/surfaces_5k/'
clubb_dir = '/home/tsw35/tyche/clubb/' # directory of clubb run folders
blank_run = '/home/tsw35/soteria/clubb/clubb_scripts/run_/' # a clean run folder to copy
cbin_dir  = '/home/tsw35/soteria/software/CLUBB/bin/' # bin for clubb
met_dir   = '/home/tsw35/soteria/clubb/data/sgp60varanarap_2012-2019.nc' #forcing
tune_o    = 'X' #'/home/tsw35/soteria/clubb/clubb_scripts/tunable_param/g0.4c110.4c80.5' 
            # 'X' overwrite tunable parameters


#### DEFAULTS ####
k       = 2 # number of clusters
nx      = 20 # number of gridcells in each direction
dx      = 5000  # resolution of surface grid in meters
dz      = 40 # dz
zmax    = 12000
dzc     = 1000
nz      = int(np.floor(zmax/dz)+1)
stdate  = '2017-07-16T11:00:00.000'# start date iso format T10:00:00.000
enddate = '2017-07-17T03:00:00.000'# end date iso format
dirname = 'test_cpl2'
c_ur    = 2
no_wind_frc = False # if true, will ignore forcings for wind velocity
l_het   = -2 # Length of heterog. in meters; -1: calculate, -2: || to wind
c_7     = -1 # default is non constant (-1), bouancy term, .3 to .8 is listed range
flux_on = 0 # turn on cirrculation flux (1) turn off (0)
wind_cr = 1 # turn on directional wind corrections to flux (1) turn off (0)
inc     = 1 # increase in vertical windspeed relative to what it should be by volume balance
cut_off = .8 # cuttoff for correlation function heterogeneity

## CIRCULATION PARAMETERS ##
# METHODS
# ri     : this uses the richardson number to determine circulation height;
#          in this case, cc1 is the richardson number cuttoff; cc1=.5
# vpt1   : this is the old vpt case and is innacurate and depricated; c11=0
# vpt2   : circulation is determined by cummulative value of gradient of vpt
#          c11 is the cuttoff value; circulation top is this cuttoff in hot
#          and circulation bottom is cuttoff in cold column.
# pa     : circulation height in hot determined by cumulative value of grad
#          of vpt. Pressure guides everything else following BJERKNES
#          circulation theorem
circ_m  = 'eb' # circulation method
cc1     = 2.8 # has different meanings depending on the circulation method
cc2     = 1 # power law for u_r
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
prs.add_argument('-r', action='store', dest='cur',type=float,default=c_ur)
prs.add_argument('-v', action='store', dest='c_7',type=float,default=c_7)
prs.add_argument('-f', action='store', dest='flux',type=int,default=flux_on)
prs.add_argument('-w', action='store', dest='wc',type=int,default=wind_cr)
prs.add_argument('-y', action='store', dest='dt',type=float,default=delta_t)
prs.add_argument('-o', action='store', dest='tun',default=tune_o)
prs.add_argument('-u', action='store', dest='cut',type=float,default=cut_off)
prs.add_argument('-a', action='store', dest='atm',default='')
prs.add_argument('-b', action='store', dest='inc',type=float,default=inc)
prs.add_argument('-cm', action='store',dest='cm',default=circ_m)
prs.add_argument('-c1', action='store',dest='cc1',type=float,default=cc1)
prs.add_argument('-c2', action='store',dest='cc2',type=float,default=cc2)

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
c_ur    = args.cur
c_7     = args.c_7
flux_on = args.flux==1
wind_cr = args.wc==1
delta_t = args.dt
tune_o  = args.tun #overwrite tunable param w/ other file
cut_off = args.cut
atm_dir = args.atm
inc     = args.inc
circ_m  = args.cm
cc1     = args.cc1
cc2     = args.cc2

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
    varout_g = np.zeros((nt,nx,nx))
    varout_v = np.zeros((nt,nx*nx))
    for t in range(nt):
        dt = stdt+datetime.timedelta(hours=t)
        tnm = dt.strftime('%Y_%m_%d_%H')
        file = sfc_dir+var+'/'+tnm+'.tif'
        data = rasterio.open(file).read(1)
        varout_g[t,:,:]=data
        varout_v[t,:]=np.reshape(data,(nx*nx))
    return varout_g,varout_v

#### CALC_BLH ####
def calc_blh(ri_,lim):
    idx=0
    try:
        wheres=np.where(ri_[:]>lim)[0]
        idx=wheres[0]
        if (idx<=3)&(len(wheres)>=2):
            idx=wheres[1]
        if (idx<=3)&(len(wheres)>=3):
            idx=wheres[2]
        if (idx<=3)&(len(wheres)>=4):
            idx=wheres[3]
    except:
        pass
    return idx

#### CLUSTER ####
def cluster(var,nt_):
    varn = normalize(var)
    ctff=0
    prev=np.percentile(varn,22.5)
    diff=100
    for i in np.linspace(25,75,21):
        ctf=np.percentile(varn,i)
        if (ctf-prev)<diff:
            diff=ctf-prev
            ctff=ctf
        prev=ctf
    Y=np.zeros((nx*nx))
    Y[varn>=ctff]=1
    k_masks=np.zeros((nt_,nx*nx))
    for t in range(nt_):
        k_masks[t,:]=Y[:]
    return k_masks

#### COMPUTE RELATIVE HUMIDITY ####
def RH(r,P,T):
    e=r*P/(.622+r)
    esat=.61078*np.exp((17.27*(T-273.15))/(T-273.15+237.3))*1000
    #print(str(e)+','+str(esat))
    return(e/esat)


#### COMPUTE CIRCULATION TENDENCIES ####
def circulate(vpt,rtm,Axx,T,H,wd,c_=c_ur,T0=T0,l=l_het,dz_=dz,dzi_=dzi,u_=[0,0],cd_=[0,0],cc1=cc1,cc2=cc2,circ_m=circ_m,pa=None,urpv=0):
    """
    from the atmospheric profiles, calculate the tendencies in moisture and
    heat as well as the vertical velocity (w) profile necessary

    Parameters
    ----------
    vpt: float(2,nz)
        Virtual Potential Temperature profiles OR Ri profiles
    rtm: float(2,nz)
        Moisture profiles
    Axx: float(2)
        vertical cross sectional area of each column
    T  : float(2,nz)
        Potential Temperature profiles
    H  : float(2)
        Area average sensible heat
    wd : float
        Width of the connection between the columns
    c_ : float
        Coefficient to reduce or increase circulation velocity
    T0 : float
        Standard Temperature
    l  : float
        Lengthscale of heterogeneity
    dz_: float
        height of each layer
    dzi: int
        depth of circulation in number of layers
    u_ : float(2)
        2 component mean wind
    cd_: float(2)
        2 component direction of circulation
    """
    if circ_m == 'eb':
        return eb_circulate(vpt,rtm,Axx,T,H,wd,c_,T0,l,dz_,dzi_,u_,cd_,cc1,cc2,circ_m,pa,urpv)

    # create empty outputs
    nz_  = vpt.shape[1]

    dThm = np.zeros((2,nz_))
    dRtm = np.zeros((2,nz_))
    wbar = np.zeros((2,nz_))

    # determine the high and low sensible heat flux column
    sgn = (H[1]-H[0])/np.abs(H[1]-H[0])
    if sgn<0:
        k_lo = 1
        k_hi  = 0
    else:
        k_lo = 0
        k_hi  = 1
    
    if circ_m == 'pa':
        # cum vpt circ height; then pressure differential
        vpthg=np.gradient(vpt[k_hi,0:int(6000/dz)])
        vpthg[0:5]=0
        vpthgc=np.cumsum(vpthg)
        maxzh=np.argmin(np.abs(vpthgc-cc1))
        pa_hm=pa[k_hi,maxzh]
        maxzl=np.argmax((pa[k_lo,:]-pa_hm)<0)-1
        minzh=np.argmax((pa[k_hi,:]-pa[k_lo,:])>0)-1
        minzl=minzh
    
        l_maxzl=min(maxzl-minzl,minzl-1)
        pa_lm=pa[k_lo,l_maxzl]
        l_maxzh=np.argmax((pa[k_hi,:]-pa_lm)<0)-1
    
    elif circ_m == 'den':
        # cum vpt circ height; then density differential
        vpthg=np.gradient(vpt[k_hi,0:int(6000/dz)])
        vpthg[0:5]=0
        vpthgc=np.cumsum(vpthg)
        maxzh=np.argmin(np.abs(vpthgc-cc1))
        vpt_hm=vpt[k_hi,maxzh]
        maxzl=np.argmax((vpt[k_lo,:]-vpt_hm)>0)-1
        try:
            minzh=np.argmax((vpt[k_hi,:]-vpt[k_lo,:])<0)-1
        except:
            print('E00: No Crossover point in VPT profiles detected')
            return dThm, dRtm, wbar, 0, [0,0], 0,0,0,0,0
        minzl=minzh

        l_maxzl=min(maxzl-minzl,minzl-1)
        l_maxzh=l_maxzl
        
    elif circ_m == 'den2':
        try:
            minzh=np.argmax((vpt[k_hi,:]-vpt[k_lo,:])<0)-1
        except:
            print('E00: No Crossover point in VPT profiles detected')
            return dThm, dRtm, wbar, 0, [0,0], 0,0,0,0,0
        minzl=minzh
        ctf_vpt=vpt[k_lo,minzh+1]-vpt[k_hi,minzh+1]
        maxzh=np.argmax((vpt[k_lo,minzh+1:-1]-vpt[k_hi,minzh+1:-1])<ctf_vpt*.5)+minzh
        maxzl=maxzh
        l_maxzl=min(maxzl-minzl,minzl-1)
        l_maxzh=l_maxzl
    
    # check positive circulation thickness for all levels
    thick = min(maxzh-minzh,maxzl-minzl,l_maxzh,l_maxzl)
    if thick<1:
        print('E01: Circulation has non-positive thickness')
        return dThm, dRtm, wbar, 0, [0,0], 0,0,0,0,0

    # create a filter that can smooth the circulation
    # adj_l is for lower, adj_u is for upper (recirc)
    a = 1.5
    b = 1.5
    adj_l = beta.pdf(np.linspace(beta.ppf(.001,a,b),beta.ppf(.999,a,b),l_maxzh),a,b)
    adj_l = adj_l/np.sum(adj_l)*(l_maxzh)

    adj_u = beta.pdf(np.linspace(beta.ppf(.001,a,b),beta.ppf(.999,a,b),maxzl-minzl),a,b)
    adj_u = adj_u/np.sum(adj_u)*np.sum(adj_l)

    adj_lch = beta.pdf(np.linspace(beta.ppf(.001,a,b),beta.ppf(.999,a,b),l_maxzl),a,b)
    adj_lch = adj_lch/np.sum(adj_lch)*np.sum(adj_l)

    adj_uhc = beta.pdf(np.linspace(beta.ppf(.001,a,b),beta.ppf(.999,a,b),maxzh-minzh),a,b)
    adj_uhc = adj_uhc/np.sum(adj_uhc)*np.sum(adj_l)

    # compute ur 
    Thi = np.mean(vpt[k_hi,0:l_maxzh])
    Tlo = np.mean(vpt[k_lo,0:l_maxzl])
    
    # check if ur is positive
    if (Thi-Tlo) < 0:
        print('E02: Surface temperature differential opposes flux diff.')
        return dThm, dRtm, wbar, 0, [0,0], 0,0,0,0,0
    
    ur = c_*((Thi-Tlo)/T0)**(cc2)*9.81**(.5)*l**(.5)
    pur = ur

    # adjust ur for mean velocity
    if not ((u_[0]==0) and (u_[1]==0)):
        circ_dir=cd_
        flow=circ_dir/np.sqrt(circ_dir[0]**2+circ_dir[1]**2)*ur
        adj_flow = flow-np.dot(np.dot(flow,u_)/np.dot(u_,u_),u_)
        ur=np.sqrt(adj_flow[0]**2+adj_flow[1]**2)

    # check if ur change is reasonable
    urpv=np.max(urpv)
    if urpv<0:
        urpv=0
    if (ur-urpv)>.5:
        ur=urpv+.5
    
    # compute horizontal circulation
    vflux = wd*dz_*ur*adj_l
    vflux_u = wd*dz_*ur*adj_u
    vflux_uhc = wd*dz_*ur*adj_uhc 
    vflux_lch = wd*dz_*ur*adj_lch

    Tl_p=np.interp(np.linspace(0,l_maxzh-1,l_maxzh),np.linspace(0,l_maxzh-1,l_maxzl),T[k_lo,0:l_maxzl])
    Rl_p=np.interp(np.linspace(0,l_maxzh-1,l_maxzh),np.linspace(0,l_maxzh-1,l_maxzl),rtm[k_lo,0:l_maxzl])
    
    dThm[k_hi,0:l_maxzh]=(vflux*Tl_p+Axx[k_hi]*dz_*T[k_hi,0:l_maxzh])/\
                      (vflux+Axx[k_hi]*dz_)-T[k_hi,0:l_maxzh]
    dRtm[k_hi,0:l_maxzh]=(vflux*Rl_p+Axx[k_hi]*dz_*rtm[k_hi,0:l_maxzh])/\
                      (vflux+Axx[k_hi]*dz_)-rtm[k_hi,0:l_maxzh]
   

    # compute recirculation
    Th_p=np.interp(np.linspace(minzl,maxzl-1,maxzl-minzl),np.linspace(minzh,maxzh-1,maxzh-minzh),T[k_hi,minzh:maxzh])
    Rh_p=np.interp(np.linspace(minzl,maxzl-1,maxzl-minzl),np.linspace(minzh,maxzh-1,maxzh-minzh),rtm[k_hi,minzh:maxzh])

    dThm[k_lo,minzl:maxzl]=dThm[k_lo,minzl:maxzl] +(vflux_u*Th_p+\
                         Axx[k_lo]*dz_*T[k_lo,minzl:maxzl])/(vflux_u+Axx[k_lo]*dz_)-\
                         T[k_lo,minzl:maxzl]
    dRtm[k_lo,minzl:maxzl]=dRtm[k_lo,minzl:maxzl] +(vflux_u*Rh_p+\
                         Axx[k_lo]*dz_*rtm[k_lo,minzl:maxzl])/(vflux_u+Axx[k_lo]*dz_)-\
                         rtm[k_lo,minzl:maxzl]
    

    # compute vertical circulation (k_hi)
    if vert_circ:
        wbar[k_hi,0:l_maxzh]=np.cumsum(vflux)/Axx[k_hi]*inc
        wbar[k_hi,l_maxzh:minzh]=wbar[k_hi,l_maxzh-1]
        wbar[k_hi,minzh:maxzh]=wbar[k_hi,l_maxzh-1]-np.cumsum(vflux_uhc)/Axx[k_hi]*inc
        wbar[k_lo,0:l_maxzl]=-np.cumsum(vflux_lch)/Axx[k_lo]*inc
        wbar[k_lo,l_maxzl:minzl]=wbar[k_lo,l_maxzl-1]
        wbar[k_lo,minzl:maxzl]=wbar[k_lo,l_maxzl-1]+np.cumsum(vflux_u)/Axx[k_lo]*inc

    urz =np.max(np.abs(wbar),axis=1)
    dzi_=(maxzh-minzh)*10**6+(maxzl-minzl)*10**3+l_maxzl+l_maxzh*10**(-3)
    print('maxzh: %f' %maxzh)
    print('maxzl: %f' %maxzl)
    print('minzh: %f' %minzh)
    print('minzl: %f' %minzl)
    print('ur: %f' %ur)
    print('urz: %f' %np.max(urz))
    


    return dThm,dRtm,wbar,ur,urz,maxzh,maxzl,dzi_,pur,Thi-Tlo

#### COMPUTE ENERGY BASED CIRCULATION ####
def eb_circulate(vpt,rtm,Axx,T,H,wd,c_,T0,l,dz_,dzi_,u_,cd_,cc1,cc2,circ_m,pa,urpv):

    # create empty outputs
    nz_  = vpt.shape[1]
    
    dThm = np.zeros((2,nz_))
    dRtm = np.zeros((2,nz_))
    wbar = np.zeros((2,nz_))
    
    k_hi=1
    k_lo=0

    sgn=(vpt[k_hi,:]-vpt[k_lo,:])/np.abs(vpt[k_hi,:]-vpt[k_lo,:]+.0000001)
    ur = c_*(np.abs(vpt[k_hi,:]-vpt[k_lo,:])/T0)**(cc2)*9.81**(.5)*l**(.5)
    pur = ur.copy()

    # adjust ur for mean velocity
    if not ((u_[0]==0) and (u_[1]==0)):
        circ_dir=np.array([sgn,sgn])
        circ_dir[0,:]=circ_dir[0,:]*cd_[0]
        circ_dir[1,:]=circ_dir[1,:]*cd_[1]
        
        flow=circ_dir/np.sqrt(circ_dir[0]**2+circ_dir[1]**2)*ur
        adj_flow = np.zeros(flow.shape)
        for x in range(nz_):
            adj_flow[:,x] = flow[:,x]-np.dot(np.dot(flow[:,x],u_)/np.dot(u_,u_),u_)
        ur=np.sqrt(adj_flow[0,:]**2+adj_flow[1,:]**2)*sgn

    # create empty energy 
    KE=np.zeros((2,len(ur)))
    KEh=np.zeros((2,len(ur)))
    BE=np.zeros((2,len(ur)))
    KEv=np.zeros((2,len(ur)))
    KE0=np.zeros((2,len(ur)))
    
    # calculate KE0
    KE0[k_hi,:]=.5*Axx[k_hi]*dz*pa[k_hi,:]**2
    KE0[k_lo,:]=.5*Axx[k_lo]*dz*pa[k_lo,:]**2

    # first level
    KEh[k_hi,0]=.5*W[0,1]*np.mean(ur[0:5])**3*dz_*s_rest
    BE[k_hi,0]=0
    KE_KEv=KEh[k_hi,0]+KE0[k_hi,0]

    for x in range(1,nz_):
        KEh[k_hi,x]=.5*W[0,1]*np.mean(ur[x])**3*dz_*s_rest
        
        # BE Calculation
        tenv=vpt[k_hi,x]/(100000/pa[k_hi,x])**.286
        rh=RH(rtm[k_hi,x-1],pa[k_hi,x-1],T[k_hi,x-1])
        if rh>=1:
            ad=9.808*(1+(2501000*r[:,x-1]/rh)/(287*T[k_hi,x-1]))/\
               (1003.5+(2501000**2*rtm[k_hi,x-1]/rh)/(461.5*T[k_hi,x-1]**2))
        else:
            ad=.0098
        tprc=(T[k_hi,x-1]-ad*dz_)*(1+.61*rtm[:,x-1])
        BE[k_hi,x]=Axx[k_hi]*9.8*(tprc-tenv)/(2*tenv)*wbar[k_hi,x-1]*dz_*s_rest
        
        # continuing...
        KE[k_hi,x]=KE0[k_hi,x]+BE[k_hi,x]+KEv[k_hi,x-1]

    '''
    # first level
    KEh[k_hi,0]=.5*W[0,1]*np.mean(ur[0:5])**3*dz_
    KEh[k_lo,0]=-.5*W[0,1]*np.mean(ur[0:5])**3*dz_
    
    KE[k_hi,0]=KEh[k_hi,0]
    KE[k_lo,0]=KEh[k_lo,0]
    wbar[k_hi,0]=np.cbrt(2*KE[k_hi,0]/Axx[k_hi])
    wbar[k_lo,0]=np.cbrt(2*KE[k_lo,0]/Axx[k_lo])

    # going up
    maxzh=0
    zcrit=0
    for x in range(1,nz_):
        if x<=5:
            KEh[k_hi,x]=.5*W[0,1]*np.mean(ur[0:5])**3*dz_
            KEh[k_lo,x]=-.5*W[0,1]*np.mean(ur[0:5])**3*dz_
        else:
            KEh[k_hi,x]=.5*W[0,1]*np.mean(ur[x])**3*dz_
            KEh[k_lo,x]=-.5*W[0,1]*np.mean(ur[x])**3*dz_
        if x>4:
            tenv=(vpt[k_hi,x-1]+vpt[k_hi,x])/2
            tprc=vpt[k_hi,x-1]
            BE[k_hi,x]=Axx[k_hi]*9.8*(tprc-tenv)/(2*tenv)*wbar[k_hi,x-1]*dz_

            tprc=vpt[k_lo,x]
            tenv=(vpt[k_lo,x-1]+vpt[k_lo,x])/2
            BE[k_lo,x]=-Axx[k_lo]*9.8*(tprc-tenv)/(2*tenv)*wbar[k_lo,x-1]*dz_
        KE[k_hi,x]=BE[k_hi,x]+KEh[k_hi,x]+KE[k_hi,x-1]
        wbar[k_hi,x]=np.cbrt(2*KE[k_hi,x]/Axx[k_hi])
        if ur[x]>0:
            KE[k_lo,x]=BE[k_lo,x]+KEh[k_lo,x]+KE[k_lo,x-1]
            if (KE[k_lo,x-1]<0) & (KE[k_lo,x]>0):
                wbar[k_lo,x]=0
                KE[k_lo,x]=0
            else:
                wbar[k_lo,x]=np.cbrt(2*KE[k_lo,x]/Axx[k_lo])
        else:
            if zcrit==0:
                zcrit=x
        

        if KE[k_hi,x]<0:
            maxzh=x-1
            KE[k_hi,x]=0
            wbar[k_hi,x]=0
            break
    ur[wbar[k_hi,:]==0]=0
    maxzl=maxzh
    
    # copy ur setting first 5 levels to positive but small
    ur_tmp=ur.copy()
    ur_tmp[0:5]=.00001

    if not zcrit==np.argmin(ur_tmp>0):
        print()
        print('zcrit error')
        print(zcrit)
        print(np.argmin(ur_tmp>0))
        print()
    zcrit=np.argmin(ur_tmp>0)

    if ur[maxzl]<0:
        # top level cold
        KEh[k_lo,maxzh]=np.abs(KEh[k_hi,maxzh])

        KE[k_lo,maxzh]=np.abs(KEh[k_lo,0])
        wbar[k_lo,maxzh]=-np.cbrt(2*KE[k_lo,0]/Axx[k_lo])
        # going down
        for x in range(maxzh-1,zcrit-1,-1):
            KEh[k_lo,x]=np.abs(KEh[k_hi,x])
            tenv=(vpt[k_lo,x+1]+vpt[k_lo,x])/2
            tprc=vpt[k_lo,x+1]
            BE[k_lo,x]=Axx[k_lo]*9.8*(tprc-tenv)/(2*tenv)*wbar[k_lo,x+1]*dz_
            KE[k_lo,x]=BE[k_lo,x]+KEh[k_lo,x]+KE[k_lo,x+1]
            if (KE[k_lo,x]<0)&(KEh[k_lo,x]>0):
                KE[k_lo,x]=0
                BE[k_lo,x]=-KEh[k_lo,x]
            wbar[k_lo,x]= -np.cbrt(2*KE[k_lo,x]/Axx[k_lo])
    
    '''
    # compute dThm and dRtm
    mskh=KEh[k_hi,:]>0
    mskl=KEh[k_lo,:]>0
    vflx=ur[:]*wd*dz_
    dThm[k_hi,mskh]=(vflx[mskh]*T[k_lo,mskh]+Axx[k_hi]*dz_*T[k_hi,mskh])/(vflx[mskh]+Axx[k_hi]*dz_)-T[k_hi,mskh]
    dThm[k_lo,mskl]=(vflx[mskl]*T[k_hi,mskl]+Axx[k_lo]*dz_*T[k_lo,mskl])/(vflx[mskl]+Axx[k_lo]*dz_)-T[k_lo,mskl]

    dRtm[k_hi,mskh]=(vflx[mskh]*rtm[k_lo,mskh]+Axx[k_hi]*dz_*rtm[k_hi,mskh])/(vflx[mskh]+Axx[k_hi]*dz_)-rtm[k_hi,mskh]
    dRtm[k_lo,mskl]=(vflx[mskl]*rtm[k_hi,mskl]+Axx[k_lo]*dz_*rtm[k_lo,mskl])/(vflx[mskl]+Axx[k_lo]*dz_)-rtm[k_lo,mskl]
    
    urz =np.max(np.abs(wbar),axis=1)

    print('maxzh: %f' %maxzh)
    print('maxzl: %f' %maxzl)
    print('ur_tp: %f' %ur[maxzh])
    print('ur_0: %f' %ur[0])
    print('ur: %f' %np.nanmax(ur))
    print('wbar: %f' %np.nanmax(wbar))

    return dThm,dRtm,wbar,ur,urz,maxzh,maxzl,dzi_,pur,0

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

# pull in the sw down and surface pressure
t0_var = datetime.datetime(2012,5,1,0,0)
fp_sw=nc.Dataset(met_dir,'r')
swt_init =int(((stdt-t0_var).total_seconds()/3600+1)) #initial time to read
swt_final=int(((endt-t0_var).total_seconds()/3600+1)) #final time to read
sw_dwn =fp_sw['sw_dn_srf'][swt_init:swt_final]
psurf  =np.nanmean(fp_sw['p_srf_aver'][swt_init:swt_final])

if np.isnan(psurf):
    psurf=97000

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
    elif line[0:5]=='p_sfc':
        psfcline='p_sfc = '+str(int(psurf))+'.e2'+'    ! Pressure at model '+\
                 'base              [Pa]'
        lines.append(psfcline+'\n')
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
        elif ((d == 'T_f[K\s]') or (d == 'rtm_f[kg\kg\s]')) and no_tq_frc:
            data2[d]=np.zeros((int(nt_frc),nz_forcing))
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

if k>1:
    k_masks=cluster(Hv[hggt,:],nt)
else:
    k_masks=np.zeros((nt,nx*nx))

# ------------------- #
# WRITE SURFACE FILES #
# ------------------- #
# transfer clustered array to surface files
print('OVERWRITE SURFACE FILES',flush=True)

# overwrite arm_sfc.in
for j in list(range(1,k+1)):
    
    sfc_path = w_dir+'/k_'+str(k)+'/c_'+str(j)+\
                   '/input/case_setups/arm_sfc.in'
    if j ==1:
        shutil.copy(sfc_path,w_dir+'/k_'+str(k)+'/arm_sfc_original.in')

    # begin overwriting file
    fp = open(sfc_path,'w')

    # Write the header information
    fp.write('! $Id$\n')
    fp.write('!\n')
    fp.write('! NoTE that T_sfc is included here,'+\
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
if k>=2:
    for i in range(nx):
        for j in range(nx):
            cl = int(clst_2d[0,i,j])
            if i>0:
                cl_2 = int(clst_2d[0,i-1,j])
                W[cl,cl_2]=W[cl,cl_2]+dx
                cdirect[cl,cl_2,1]=cdirect[cl,cl_2,1]-1
            if i<(nx-1):
                cl_2 = int(clst_2d[0,i+1,j])
                W[cl,cl_2]=W[cl,cl_2]+dx
                cdirect[cl,cl_2,1]=cdirect[cl,cl_2,1]+1
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

# compute the area of each cluster
A = clst_frac*nx*nx*dx*dx
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
min_rad = 400
d_frac = s_rest/decay_t # percent change per circ timestep in initiation/decay
fluxbool=np.zeros((len(sw_dwn_interp),))
fluxbool[sw_dwn_interp>min_rad]=1
if flux_on:
    # change boolean to decay setup
    doflux=[0]
    for i in range(1,len(fluxbool)):
        prev=doflux[i-1]
        if (prev<.99)&(fluxbool[i]==1):
            doflux.append(prev+d_frac)
        elif (prev>.01)&(fluxbool[i]==0):
            doflux.append(prev-d_frac)
        elif fluxbool[i]==0:
            doflux.append(0)
        elif fluxbool[i]==1:
            doflux.append(1)
        else:
            print('ERROR -- FLUXBOOL NOT 0 or 1')
else:
    doflux=np.zeros((len(sw_dwn_interp),))
print()
print('DOFLUX')
print(doflux)
print()

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
fp.write('T0      = '+str(T0)+' # clubb reference temperature\n')
fp.write('c_ur    = '+str(c_ur)+' # circulation velocity coefficient\n')
fp.write('wind frc= '+str(no_wind_frc==False)+' # includes wind forcings?\n')
fp.write('l_het   = '+str(l_het)+' # lengthscale of heterogeneity\n')
fp.write('sfc_dir = '+sfc_dir+' # sfc directory\n')
fp.write('k       = '+str(k)+' # number of columns\n')
fp.write('nx      = '+str(nx)+' # number of gridcells in each direction\n')
fp.write('dx      = '+str(dx)+' # surface grid resolution\n')
fp.write('dz      = '+str(dz)+' # vertical grid resolution\n')
fp.write('zmax    = '+str(zmax)+' # maximum height in clubb\n')
fp.write('dzc     = '+str(dzc)+' # thickness of circulation\n')
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
dvpt = np.ones((len(tlist),))*-1
u_r = np.ones((len(tlist),nz))*-1
u_r0 = np.ones((len(tlist),nz))*-1
u_rz=np.ones((2,len(tlist)))*-1
z_circh = np.ones((len(tlist),))*-1
z_circl = np.ones((len(tlist),))*-1
dzic = np.ones((len(tlist),))*-1
wz_s = np.zeros((k,len(tlist),nz))


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
    rig  = np.zeros((k,nz))
    pa   = np.zeros((k,nz))
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
            if circ_m== 'ri':
                fp2 = nc.Dataset(m_dir+'/c_'+str(j)+'/output/arm_zm.nc','r')
                rig[j-1,:]=fp2['Richardson_num'][int(round(n_rest_out-1)),:,0,0]
                fp2.close()
            if circ_m== 'pa':
                pa[j-1,:]=fp['p_in_Pa'][int(round(n_rest_out-1)),:,0,0]
            if circ_m== 'eb':
                pa[j-1,:]=fp['p_in_Pa'][int(round(n_rest_out-1)),:,0,0]

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
        dThlm=np.zeros((k,thvm.shape[1]))
        dRtm=np.zeros((k,thvm.shape[1]))
        for k1 in range(k):
            for k2 in range(k1+1,k):
                if circ_m=='eb':
                    urpvin=wz_s[:,i-1,:]
                else:
                    urpvin=u_r[i-1,:]
                if circ_m=='ri':
                    data=rig[[k1,k2],:]
                else:
                    data=thvm[[k1,k2],:]
                dThlm,dRtm,wz_s[:,i,:],u_r[i,:],u_rz[:,i],z_circh[i],z_circl[i],dzic[i],u_r0[i,:],dvpt[i] = \
                        circulate(data,rtms[[k1,k2],:],A[[k1,k2]],\
                        tmps[[k1,k2],:],H2[i,[k1,k2]],W[k1,k2],c_=c_ur*doflux[i],\
                        cd_=cdirect[k1,k2,:],u_=um,l=l_het,pa=pa,urpv=urpvin)

    for j in list(range(1,k+1)):
        frc_file = m_dir+'/c_'+str(j)+'/input/case_setups/arm_forcings.in'
       
        # Load Tendencies based on fluxes
        if doflux[i]>.01:
            #frcs[j]['T_f[K\s]'][i,:]      = frcs[j]['T_f[K\s]'][i,:] + dThlm[j-1,:]
            #frcs[j]['rtm_f[kg\kg\s]'][i,:]=frcs[j]['rtm_f[kg\kg\s]'][i,:] + dRtm[j-1,:]
            frcs[j]['T_f[K\s]'][i+1,:]      = frcs[j]['T_f[K\s]'][i+1,:] + dThlm[j-1,:]
            frcs[j]['rtm_f[kg\kg\s]'][i+1,:]=frcs[j]['rtm_f[kg\kg\s]'][i+1,:] + dRtm[j-1,:]

            # change from m/s to Pa/s
            omega=w_to_omega(wz_s[j-1,i,:],frcs[j]['Press[Pa]'][i,:],tmps[j-1,:])
            #frcs[j]['omega[Pa\s]'][i,:]=frcs[j]['omega[Pa\s]'][i,:]+omega
            frcs[j]['omega[Pa\s]'][i+1,:]=frcs[j]['omega[Pa\s]'][i+1,:]+omega

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
        cmd = rfile+' arm -p '+cbasedir+'/input/tunable_parameters.in'
        
        print('STARTING RUN: column '+str(j)+' at time '+str(t0/60/60-5)+':'+str(t0),flush=True)
        try:
            subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL)
        except Exception as e:
            print('ERROR RUN FAILED')
            print(e)
            print()
            break
        print('RUN COMPLETE: column '+str(j),flush=True)
    print()

#### OUTPUT ####
# Output Final forcings
for j in list(range(1,k+1)):
    write_forcings(frcs[j],w_dir+'/k_'+str(k)+'/c_'+str(j)+'/arm_forcings_f.in')

# output circulation characteristics
fp=nc.Dataset(w_dir+'/k_'+str(k)+'/clusters.nc','r+')
fp.createDimension('z',size=nz)
fp.createDimension('xy',size=2)
fp.createDimension('t_model',size=len(tlist))
fp.createVariable('t_model','d',dimensions=('t_model'))
fp.createVariable('u_r','d',dimensions=('t_model','z'))
fp.createVariable('z_circh','d',dimensions=('t_model'))
fp.createVariable('z_circl','d',dimensions=('t_model'))
fp.createVariable('dz_circ','d',dimensions=('t_model'))
fp.createVariable('flow_dir','d',dimensions=('clst','clst','xy'))
fp.createVariable('u_rz','d',dimensions=('clst','t_model'))
fp.createVariable('u_r0','d',dimensions=('t_model','z'))
fp.createVariable('dvpt','d',dimensions=('t_model'))
fp.createVariable('wbar','d',dimensions=('clst','t_model','z'))
fp['u_r'][:]=u_r[:]
fp['z_circh'][:]=z_circh[:]*dz
fp['z_circl'][:]=z_circl[:]*dz
fp['dz_circ'][:]=dzic[:]
fp['flow_dir'][:]=cdirect[:]
fp['u_rz'][:]=u_rz[:]
fp['u_r0'][:]=u_r0[:]
fp['dvpt'][:]=dvpt[:]
fp['wbar'][:]=wz_s[:]

        

    
    

    

