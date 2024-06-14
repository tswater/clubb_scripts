# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import datetime
import scipy as sci
import scipy.ndimage
import re
import os
import rasterio
import seaborn as sns
import matplotlib as mpl
from scipy.stats import beta
from sklearn.linear_model import LinearRegression
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import ImageGrid
mpl.rcParams['figure.dpi'] = 300
sns.set_theme()
plt.rcParams.update({'figure.max_open_warning': 0})



# %%
### DIRECTORIES ###
clubb_dir='/home/tsw35/tyche/clubb/'
H0dir=clubb_dir+'sgp_cpl_l1_w/'
subdir= 'may111_dz60dt06'#'apr25_dz60dt06'
subdir= 'may24_dz60dt06'
les_1c   ='/home/tsw35/tyche/data/LES_1C/'
les_prmd ='/home/tsw35/soteria/clubb/data/les_param2/'
les_prmd1='/home/tsw35/soteria/clubb/data/les_param/'
tskdir='/home/tsw35/soteria/clubb/data/surfaces_5k/lw'
sfc_dir   = '/home/tsw35/soteria/clubb/data/surfaces_5k/'
dir_cpl  =clubb_dir+subdir+'_cpl/'
dir_1c   =clubb_dir+subdir+'_1c/'
dir_2c   =clubb_dir+subdir+'_2c/'
prime_day='20170802'
core_days=['20160625','20170717','20180709'] #['20160625','20170716','20170717']#['20160625','20170717','20180709']
dayst={}
for day in core_days:
    dayst[day]=day[0:4]+'-'+day[4:6]+'-'+day[6:8]
km5=83
km5l=166    


# %% [markdown]
# # HELPER FUNCTIONS 

# %%
def read_sfc_data(var,nt,stdt,override='X',nx=20):
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

def circ_h(vpt_,dtsk_,cc1=.85):
    maxzh=np.zeros((16,))
    minzh=np.zeros((16,))
    maxzl=np.zeros((16,))
    for t in range(1,16):
        vpthg=np.gradient(vpt_[t,1,:])
        vpthg[0:5]=0
        cdT=cc1*dtsk_[t,0]
        vpthgc=np.cumsum(vpthg)
        maxzh[t]=np.argmin(np.abs(vpthgc-cdT))
        vpt_hm=vpt[t,1,int(maxzh[t])]
        maxzl[t]=np.argmax((vpt_[t,0,:]-vpt_hm)>0)-1
        try:
            minzh[t]=np.where((vpt_[t,1,:]-vpt_[t,0,:])<0)[0][0]-1
        except:
            minzh[t]=0
    return maxzh,maxzl,minzh

def circ2thick(cc,mzz,test=False):
    cc2=cc.copy()
    for t in range(cc.shape[1]):
        mz=mzz[t].astype(int)
        #low=np.where(np.isnan(cc[:,t]))[0][0]
        low=np.where(cc[:,t]==0)[0][0]
        if low==0:
            continue
        #mn=np.where(~np.isnan(cc[:,t]))[0][low]
        mn=np.where(cc[:,t]>0)[0][low]
        #mx=np.where(~np.isnan(cc[:,t]))[0][-1]+1
        mx=np.where(cc[:,t]>0)[0][-1]+1
        rng1=mx-mn
        rng2=mz-mn
        if rng1<=0:
            continue
        if rng2<=0:
            continue
        a = 1.5
        b = 1.5
        adj = beta.pdf(np.linspace(beta.ppf(.001,a,b),beta.ppf(.999,a,b),rng2),a,b)
        adj=adj/np.sum(adj)*np.sum(cc[mn:mx,t])
        
        cc2[mn:mz,t]=adj[:]
    return cc2

def read_forcings(path,nz_=301):
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


# %% [markdown]
# # LOAD FILES

# %%
### LOAD IN FILES ###
fcs  = {}
fss  = {}
fts  = {}
fms  = {}

fs1s  = {}
ft1s  = {}
fm1s  = {}

fs2s  = {}
ft2s  = {}
fm2s  = {}

fs2hs = {}
fs2cs = {}

fshs  = {}
fths  = {}
fmhs  = {}

fscs  = {}
ftcs  = {}
fmcs  = {}

flts = {}
flms = {}
fps  = {}
fps2  = {}
fv   = nc.Dataset('/home/tsw35/soteria/clubb/data/sgp60varanarap_2012-2019.nc','r')
for day in core_days:
    fss[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/agg_outsfc.nc','r')
    fts[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/agg_outzt.nc','r')
    fms[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/agg_outzm.nc','r')
    
    fs1s[day]=nc.Dataset(dir_1c+'sgp_'+day+'/k_1/c_1/output/arm_sfc.nc','r')
    ft1s[day]=nc.Dataset(dir_1c+'sgp_'+day+'/k_1/c_1/output/arm_zt.nc','r')
    fm1s[day]=nc.Dataset(dir_1c+'sgp_'+day+'/k_1/c_1/output/arm_zm.nc','r')
    
    fshs[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/c_2/output/arm_sfc.nc','r')
    fths[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/c_2/output/arm_zt.nc','r')
    fmhs[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/c_2/output/arm_zm.nc','r')
    
    fscs[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/c_1/output/arm_sfc.nc','r')
    ftcs[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/c_1/output/arm_zt.nc','r')
    fmcs[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/c_1/output/arm_zm.nc','r')
    
    fs2s[day]=nc.Dataset(dir_2c+'sgp_'+day+'/k_2/agg_outsfc.nc','r')
    ft2s[day]=nc.Dataset(dir_2c+'sgp_'+day+'/k_2/agg_outzt.nc','r')
    fm2s[day]=nc.Dataset(dir_2c+'sgp_'+day+'/k_2/agg_outzm.nc','r')
    
    fs2cs[day]=nc.Dataset(dir_2c+'sgp_'+day+'/k_2/c_1/output/arm_sfc.nc','r')
    fs2hs[day]=nc.Dataset(dir_2c+'sgp_'+day+'/k_2/c_2/output/arm_sfc.nc','r')
    
    flts[day]=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    flms[day]=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    
    fcs[day]=nc.Dataset(dir_cpl+'sgp_'+day+'/k_2/clusters.nc','r')
    fps[day]=nc.Dataset(les_prmd+day+'.nc','r')
    fps2[day]=nc.Dataset(les_prmd1+day+'.nc','r')

# %%
# forcings 
frc0={}
frcc={}
frch={}
for day in core_days:
    frc0[day]=read_forcings(dir_2c+'sgp_'+day+'/k_2/c_2/arm_forcings_f.in',201)
    frcc[day]=read_forcings(dir_cpl+'sgp_'+day+'/k_2/c_1/arm_forcings_f.in',201)
    frch[day]=read_forcings(dir_cpl+'sgp_'+day+'/k_2/c_2/arm_forcings_f.in',201)

# %%
frc0[day].keys()

# %%

# %% [markdown]
# # 

# %% [markdown]
# # PROFILES for DAYS

# %%
#### PROFILES for DAYS ####
#### PROFILES for DAYS ####
#### PROFILES for DAYS ####

# %%
#### PROFILES for DAYS ####
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
def index_time(date,clb_st=7,les_st=7,var_st=datetime.datetime(2012,5,1,0,0)):
    clb_i=(date.hour-clb_st)*60+date.minute
    if clb_i<0:
        clb_i=clb_i+60*24
    les_i=int((date.hour-les_st)*6+date.minute/10)
    var_i=int((date-var_st).total_seconds()/60/60)+5
    return clb_i,les_i,var_i
def vpt(T,r,p):
    theta=T*(100000/(p))**(2/7)
    thv=theta#*(1+.61*r)
    return thv


# %%
#### PROFILES for DAYS ####
# Temperature
# Moisture
# Surface
# Wind
wr=[1.75,.5,1,1,1]
fig,ax=plt.subplots(3,5,figsize=(8,10),gridspec_kw={'width_ratios':wr})
ax=ax.flatten()
#fig.tight_layout()
i=0
vmin=50
vmax=475
labels=['a)','','b)','c)','d)','e)','','f)','g)','h)','i)','','j)','k)','l)']
for day in core_days:
    top=10000
    
    altl=flts[day]['z'][1,:]
    altc=ft1s[day]['altitude'][:]
    
    st_time= datetime.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),7)
    clb_t,les_t,var_t=index_time(st_time)
    
    var_p=fv['lev'][:]
    altv=np.flip(np.interp(np.flip(var_p*100),np.flip(ft2s[day]['p_in_Pa'][clb_t,:,0,0]),np.flip(altc)))
    
    u  = flms[day]['u'][1+5*6,:]
    v  = flms[day]['v'][1+5*6,:]
    
    ax[i+1].remove()
    
    pn=3
    T  = flms[day]['thv'][1,:]
    ax[i+pn].plot(T,altl/1000,'-',c='pink',ms=3)
    T  = flms[day]['thv'][1+6*5,:]
    ax[i+pn].plot(T,altl/1000,'-',c='indianred',ms=3)
    T  = flms[day]['thv'][1+6*10,:]
    ax[i+pn].plot(T,altl/1000,'-',c='firebrick',ms=3)
    
    ax[i+pn].set_xlim([299,339])
    ax[i+pn].set_ylim([-.25,top/1000])
    ax[i+pn].set_yticks([0,2,4,6,8,10])
    ax[i+pn].set_yticklabels([])
    ax[i+pn].set_xticks([300,310,320,330])
    
    if i<5:
        ax[i+pn].legend(['7:00','12:00','17:00'],handlelength=0.5)
    
    if i >=10:
        ax[i+pn].set_xticklabels([300,310,320,330],rotation=45)
        ax[i+pn].set_xlabel(r'$\theta$ ($K$)')
    else:
        ax[i+pn].set_xticklabels([])
    ax[i+pn].set_title(labels[i+pn],loc='left')
        
        
    pn=4
    
    q  = flms[day]['qv'][1,:]
    ax[i+pn].plot(q*1000,altl/1000,'-',c='cornflowerblue',ms=3)
    q  = flms[day]['qv'][1+6*5,:]
    ax[i+pn].plot(q*1000,altl/1000,'-',c='royalblue',ms=3)
    q  = flms[day]['qv'][1+6*10,:]
    ax[i+pn].plot(q*1000,altl/1000,'-',c='mediumblue',ms=3)
    
    ax[i+pn].set_xlim([-1,18])
    if i<5:
        ax[i+pn].legend(['7:00','12:00','17:00'],handlelength=0.5)
    ax[i+pn].set_yticks([0,2,4,6,8,10])
    ax[i+pn].set_yticklabels([])
    ax[i+pn].set_xticks([0,5,10,15])
    ax[i+pn].set_title(labels[i+pn],loc='left')
    ax[i+pn].set_ylim([-.25,top/1000])
    if i >=10:
        ax[i+pn].set_xticklabels([0,5,10,15],rotation=45)
        ax[i+pn].set_xlabel('$q_t$ ($g\ kg^{-1}$)')
    else:
        ax[i+pn].set_xticklabels([])
    
    pn=2
    ax[i+pn].plot([0,0],[-1,15],'-',lw=3,c='white',label=None)
    ax[i+pn].plot(u,altl/1000,'-',c='black',ms=3,label='$u$')
    ax[i+pn].plot(v,altl/1000,':',c='black',ms=3,label='$v$')
    ax[i+pn].plot([-15,15],[0,0],'-',lw=3,c='white',label=None)
    ax[i+pn].set_xlim([-14.5,14.5])
    ax[i+pn].set_ylim([-.25,top/1000])
    ax[i+pn].set_yticks([0,2,4,6,8,10])
    ax[i+pn].set_xticks([-10,-5,0,5,10])
    ax[i+pn].set_ylabel('Altitude ($km$)')
    ax[i+pn].set_title(labels[i+pn],loc='left')
    if i<5:
        ax[i+pn].legend(handlelength=1)
    if i >=10:
        ax[i+pn].set_xticklabels([-10,-5,0,5,10],rotation=45)
        ax[i+pn].set_xlabel('Wind ($m\ s^{-1}$)')
    else:
        ax[i+pn].set_xticklabels([])
    
    im=ax[i].imshow(fcs[day]['H'][5,:,:],cmap='coolwarm',vmin=vmin,vmax=vmax)
    ax[i].grid(False)
    ax[i].set_yticklabels([])
    ax[i].set_xticklabels([])
    ax[i].set_ylabel(dayst[day],fontsize=15)
    ax[i].set_title(labels[i],loc='left')
    if i >=10:
        ax[i].set_xlabel('\n\nSensible Heat ($W\ m^{-2}$)')
    axins = inset_axes(ax[i],
                    width="100%",  
                    height="5%",
                    loc='lower center',
                    borderpad=-1)
    plt.colorbar(im,cax=axins,orientation='horizontal')
    
    '''
    im=ax[i+1].imshow(fcs[day]['LE'][7,:,:],cmap='summer')
    ax[i+1].grid(False)
    ax[i+1].set_yticklabels([])
    ax[i+1].set_xticklabels([])
    axins = inset_axes(ax[i+1],
                    width="100%",  
                    height="5%",
                    loc='lower center',
                    borderpad=-1)
    plt.colorbar(im,cax=axins,orientation='horizontal')
    '''
    i=i+5
    
plt.subplots_adjust(hspace=.15)
plt.show()

# %% [markdown]
# # 

# %% [markdown]
# # FAKE PROFILE for DIAGRAM

# %%
'''#### FAKE PROFILE for DIAGRAM ####
   #### FAKE PROFILE for DIAGRAM ####
   #### FAKE PROFILE for DIAGRAM ####'''

# %%
#### FAKE PROFILE for DIAGRAM ####
# Make fake data
ztop =800
zmax1=720
zmax2=640
zcrit=475
zcirc=240
alt=np.linspace(0,800,80)
vptw=np.zeros((80,))
vptc=np.zeros((80,))
vptc[0:32]=1
vptw[0:40]=2
vptw[36:48]=12*(.025*np.linspace(0,11,12))**2+2
vptc[28:48]=2.4*(.05*np.linspace(0,19,20))**4+1
vptw[48]=3
vptw[49:72]=10*(.025*np.linspace(10,24,23))**2+2.5
vptc[48:50]=2.4*(.05*np.linspace(20,21,2))**4+1
#vptw[72:80]=10*(.025*np.linspace(24,30,8))**2+2.5
vptw[72:74]=10*(.025*np.linspace(25,26,2))**2+2.5

#vptw[39:72]=6*(.025*np.linspace(0,32,33))**2+2
vptw[74:80]=.98*(np.sqrt(np.linspace(7.5+4,14.2,6)))+vptc[49]*.9
#vptc[31:51]=1.5*(.05*np.linspace(0,19,20))**4+1
vptc[50:64]=.98*(np.sqrt(np.linspace(.5,7,14)))+vptc[49]*.9
vptc[64:80]=.98*(np.sqrt(np.linspace(7.5,14,16)))+vptc[49]*.9

# %%
#### FAKE PROFILE for DIAGRAM ####
# Make fake plot
ylabels=[r'$0$',r'$z_{circ}$',r'$z_{crit}$',r'$z_{lim_2}$',r'$z_{lim_1}$']
xlabels=[r'$\theta_{sfc_2}$',r'$\theta_{sfc_1}$',r'$\theta_{crit}$',r'$\theta_{lim}$']
plt.figure(figsize=(3,5),dpi=300)
plt.plot(vptc,alt,c='cornflowerblue')
plt.plot(vptw,alt,c='indianred')
plt.xlabel('Virtual Potential Temperature (K)')
plt.yticks([0,zcirc,zcrit,zmax2,zmax1],labels=ylabels)
plt.xticks([1,2,2.9,6.1],labels=xlabels)
plt.show()

# %% [markdown]
# # 

# %% [markdown]
# # LES CIRCULATION MATCHING STRUCTURE/VELOCITY

# %%
#### LES CIRCULATION MATCHING ####
#### LES CIRCULATION MATCHING ####
#### LES CIRCULATION MATCHING ####

# %%
#### LES CIRCULATION MATCHING VELOCITY ####
# Velocity Matching Prep.
filelist=os.listdir(les_prmd)
filelist.sort()
u2s=[]
ums=[]
dvpts=[]
for file in filelist:
    if '20180811' in file:
        continue
    if '20160819' in file:
        continue
    if '20170705' in file:
        continue
    if '201908' in file:
        continue
    t1=3
    t2=12
    ht=16 #first 500m is 17
    ht0=5
    
    fp=nc.Dataset(les_prmd+file,'r')
    
    day=file[0:8]
    
    u90=fp['u90'][t1:t2,ht0:ht]
    u10=fp['u10'][t1:t2,ht0:ht]
    v90=fp['v90'][t1:t2,ht0:ht]
    v10=fp['v10'][t1:t2,ht0:ht]
    lhet=fp['lhet'][0]
    if lhet<0:
        print(file)
    dvpt=fp['vpt'][t1:t2,1,ht0:ht]-fp['vpt'][t1:t2,0,ht0:ht]
    
    msk=((u90/np.abs(u90))==(u10/np.abs(u10)))|(dvpt<0)
    u90[msk]=0
    u10[msk]=0
    meanu=np.abs((u90+u10)/2)
    m2=(np.abs(u10)>np.abs(u90))
    u90[m2]=u10[m2]
    u90=np.abs(u90)
    
    msk=((v90/np.abs(v90))==(v10/np.abs(v10)))|(dvpt<0)
    v90[msk]=0
    v10[msk]=0
    meanv=np.abs((v90+v10)/2)
    m2=(np.abs(v10)>np.abs(v90))
    v90[m2]=v10[m2]
    v90=np.abs(v90)
    
    u2=np.sqrt(v90**2+u90**2)
    dvpt=dvpt[u2>0]
    u2=u2[u2>0]
    #u2=u2[u2>0]*300/(lhet**(1/2)*9.8**(1/2))
    u2s.extend(u2)
    um=dvpt/300*(lhet**(1/2)*9.8**(1/2))
    ums.extend(um)
    dvpts.extend(dvpt)
        
dvpts=np.array(dvpts)
ums=np.array(ums)
u2s=np.array(u2s)

model=LinearRegression(fit_intercept=False)
msk=(ums<=20)&(u2s<=20)
model.fit(ums[msk].reshape(-1, 1),u2s[msk].reshape(-1, 1))
print(model.coef_)
print(model.score(ums[msk].reshape(-1, 1),u2s[msk].reshape(-1, 1)))
print(model.intercept_)

# %%
# Structure Matching data prep
tsk={}
dx=5000
for day in core_days:
    file=day+'.nc'
    stdt=datetime.datetime(int(file[0:4]),int(file[4:6]),int(file[6:8]),12,0)
    lwg,lwv=read_sfc_data('lw',16,stdt)
    tskg=(lwg/(5.67*10**(-8)))**(1/4)
    tsk[day]=np.array([np.std(tskg,axis=(1,2)),np.std(tskg,axis=(1,2))])

# %%
cc1=.85
frc=1.5
km4l=133
fig=plt.figure(figsize=(12,7),dpi=300)
gs = GridSpec(6, 5, figure=fig)
i=0
labels=['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)']

for day in core_days:
    file=day+'.nc'
    u90=fps[day]['u90'][:]
    u10=fps[day]['u10'][:]
    v90=fps[day]['v90'][:]
    v10=fps[day]['v10'][:]
    
    vpt=fps[day]['vpt'][:]
    dtsk=2*tsk[day].T
    
    ax1 = fig.add_subplot(gs[i*2:i*2+2, 0:3])
    
    maxzh,maxzl,minzh = circ_h(vpt,dtsk,cc1)

    data=np.sqrt(((u90[:,0:km4l].T-u10[:,0:km4l].T)/2)**2+((v90[:,0:km4l].T-v10[:,0:km4l].T)/2)**2)
    
    data=data/np.max(data,axis=0)
    
    im=ax1.imshow(data,cmap='terrain',origin='lower',extent=(7,22,0,frc*4))
    fig.colorbar(im,ax=ax1,label=r'Normalized $u_r$')
    
    ax1.plot(np.linspace(7,22,16),maxzh*30*frc/1000,'k--')
    ax1.plot(np.linspace(7,22,16),maxzl*30*frc/1000,'k-')
    ax1.plot(np.linspace(7,22,16),minzh*30*frc/1000,'k:')
    ax1.set_ylim(0,frc*4)
    ax1.set_xlim(8,20)
    ax1.set_title(labels[i],loc='left')
    ax1.set_yticks([frc,frc*2,frc*3,frc*4])
    ax1.set_yticklabels([1,2,3,4])
    ax1.set_xticks([8,12,16,20])
    if i==2:
        ax1.set_xticklabels(['8:00','12:00','16:00','20:00'])
        ax1.set_xlabel('Local Time')
    else:
        ax1.set_xticklabels([])
    if i == 1:
        ax1.set_ylabel('Altitude ($km$)',fontsize=16)
    if i == 0:
        ax1.legend([r'$z_{lim_1}$',r'$z_{lim_2}$',r'$z_{crit}$'],handlelength=1,handletextpad=.6,loc='upper left',framealpha=.9)
    
    i=i+1
ax=fig.add_subplot(gs[1:-1,3:5])
ax.plot(np.linspace(0,7),np.linspace(0,7),'--',color='grey',zorder=0)
ax.scatter(ums[msk]*1.36,u2s[msk],s=6,alpha=.1,zorder=1)
ax.set_xlim(0,7)
ax.set_title(labels[i],loc='left')
ax.set_ylim(0,7)
ax.set_xlabel('Model Velocity ($ms^{-1}$)',fontsize=16)
ax.set_ylabel('LES Velocity ($ms^{-1}$)',fontsize=16)
fig.subplots_adjust(wspace=.35,hspace=.4)
plt.show()

# %% [markdown]
# # 

# %% [markdown]
# # CIRCULATION COMPARISON LES and CLUBB

# %%
#### CIRCULATION COMPARISON LES and CLUBB ####
#### CIRCULATION COMPARISON LES and CLUBB ####
#### CIRCULATION COMPARISON LES and CLUBB ####
km5=83
km5l=166
frc=1.25
cc1=.85
vmax=4
vmin=0
testdir=clubb_dir+subdir+'_cpl/'
cmap='terrain'
i=1

plt.figure(figsize=(10,4),dpi=300)

#filelist=['20160625.nc','20170716.nc','20170717.nc','20190707.nc','20180709.nc','20180707.nc','20170802.nc','20150801.nc','20160610.nc','20160716.nc']
for day in core_days:
    file='sgp_'+day
        
    day=file[4:]
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    ft=nc.Dataset(testdir+file+'/k_2/agg_outzt.nc','r')
    
    fp=nc.Dataset(les_prmd+day+'.nc','r')

    runpm=open(testdir+file+'/tw_run_param.txt','r').readlines()
    
    for ln in runpm:
        if 'c_ur' in ln:
            c_ur=re.findall('\d+\.\d+',ln)[0]
            
    ur0 =fc['u_r0'][:,:km5]
    uu  =ft['um'][:,:km5,0,0]
    vv  =ft['vm'][:,:km5,0,0]
    wfrac=fc['W_frac'][:]
    ur=fc['u_r'][:,:km5]
    ur[ur==-1]=0#float('nan')
    ur[ur==0]=0#float('nan')
    ur=ur[0:-12,:]
    maxzh=fc['z_circh'][:-12]
    data= np.abs(ur.T)
    if day=='20160625':
        data=data/fc['W_frac'][0]
    
    data=circ2thick(data,maxzh/60,day=='20170717')
    
    plt.subplot(2,3,i)
    #data[data==0]=float('nan')
    plt.imshow(data,cmap=cmap,origin='lower',vmin=vmin,vmax=vmax,extent=(7,22,0,5*frc))
    #if i==3:
    #    plt.colorbar(label=r'$u_r$ ($ms^{-1}$)')
    #plt.title(dayst[day])
    plt.title(dayst[day]+'\n'+labels[i-1]+'                               ',fontsize=15)
    
    #plt.plot(np.linspace(7,22,len(maxzh)),maxzh/1000*frc,'k--')
    plt.ylim(0,frc*5)
    plt.xlim(8,20)
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([9,12,15,18],labels=[])
    if i==1:
        plt.ylabel('Altitude ($km$)                        ',fontsize=15)
    
    #plt.xlabel('Local Time')
    
    
    
    
    
    plt.subplot(2,3,i+3)
    u90=fps[day]['u90'][:]
    u10=fps[day]['u10'][:]
    v90=fps[day]['v90'][:]
    v10=fps[day]['v10'][:]
    
    vpt=fps[day]['vpt'][:]
    dtsk=2*tsk[day].T
    
    
    maxzh,maxzl,minzh = circ_h(vpt,dtsk,cc1)
    
    #plt.subplot(1,2,1)
    #plt.title(dayst[day])
    #data=u90[:,0:166].T-u10[:,0:166].T
    data=np.sqrt(((u90[:,0:166].T-u10[:,0:166].T)/2)**2+((v90[:,0:166].T-v10[:,0:166].T)/2)**2)
    
    nch1=(u90[:,0:km5l].T/np.abs(u90[:,0:km5l].T))==(u10[:,0:km5l].T/np.abs(u10[:,0:km5l].T))
    nch2=(v90[:,0:km5l].T/np.abs(v90[:,0:km5l].T))==(v10[:,0:km5l].T/np.abs(v10[:,0:km5l].T))
    
    
    '''
    data[nch1&nch2]=float('nan')
    
    if day=='20160625':
        data[int(3.5/5*km5l):,:]=float('nan')
        
    if day=='20170716':
        data[int(3/5*km5l):,:]=float('nan')
        
    if day=='20170717':
        data[int(3.1/5*km5l):,:]=float('nan')
    '''
    
    #data=data/np.max(data,axis=0)
    #data[data<1]=float('nan')
    plt.imshow(data,cmap=cmap,origin='lower',vmin=vmin,vmax=vmax,extent=(7,22,0,5*frc))
    #plt.colorbar(label=r'$u_r$ ($ms^{-1}$)')
    if i==3:
        plt.colorbar(label=r'$u_r$ ($ms^{-1}$)')
    
    #plt.plot(np.linspace(7,22,16),maxzh*30*frc/1000,'k--')
    #plt.plot(np.linspace(7,22,16),maxzl*30*frc/1000,'k-')
    #plt.plot(np.linspace(7,22,16),minzh*30*frc/1000,'k:')
    plt.ylim(0,frc*5)
    plt.xlim(8,20)
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([9,12,15,18],labels=['9:00','12:00','15:00','18:00'])
    if i==2:
        plt.xlabel('Local Time',fontsize=15)
    #plt.legend([r'$z_{max_1}$'],loc='upper left',framealpha=.9)
    
    i=i+1
plt.show()    

# %% [markdown]
# # 

# %% [markdown]
# # MEAN PROFILES COMPARISON

# %%
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####
plt.figure(figsize=(6,8))
alt=ft1s[day]['altitude'][0:km5+10]/1000
i=1
t_i=60*(17-7)
for day in core_days:
    plt.subplot(2,3,i)
    t1=ft1s[day]['thlm'][t_i,0:km5+10,0,0]
    t2=ft2s[day]['thlm'][t_i,0:km5+10,0,0]
    tc=fts[day]['thlm'][t_i,0:km5+10,0,0]
    plt.plot(t1,alt,':',c='indianred')
    plt.plot(t2,alt,'--',c='indianred')
    plt.plot(tc,alt,c='indianred')
    
    plt.title(labels[i-1],fontsize=14)
    plt.xlabel(r'$\theta$ ($K$)',fontsize=14)
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['SC','IC','TCM'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    plt.xlim(306,326)
    plt.xticks([310,315,320,325],rotation=30)
    plt.ylim(0,5.1)
    
    
    plt.subplot(2,3,i+3)
    q1=ft1s[day]['rtm'][t_i,0:km5+10,0,0]*1000
    q2=ft2s[day]['rtm'][t_i,0:km5+10,0,0]*1000
    qc=fts[day]['rtm'][t_i,0:km5+10,0,0]*1000
    plt.plot(q1,alt,':',c='cornflowerblue')
    plt.plot(q2,alt,'--',c='cornflowerblue')
    plt.plot(qc,alt,c='cornflowerblue')
    
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['SC','IC','TCM'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    plt.xlabel('$q_t$ ($g\ kg^{-1}$)',fontsize=14)
    plt.xticks([0,5,10],['','5','10'])
    plt.xlim(-1,13)
    plt.ylim(0,5.1)
    
    i=i+1
    
plt.subplots_adjust(hspace=.3)
plt.show()

# %%
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####
plt.figure(figsize=(6,8))
alt=flts[day]['z'][50,0:km5l+10]/1000
i=1
t_i=6*(17-7)
for day in core_days:
    plt.subplot(2,3,i)
    t1=flms[day]['thl'][t_i,0:km5l+10]
    tc=flts[day]['thl'][t_i,0:km5l+10]
    plt.plot(t1,alt,':',c='indianred')
    plt.plot(tc,alt,c='indianred')
    
    plt.title(labels[i-1],fontsize=14)
    plt.xlabel(r'$\theta$ ($K$)',fontsize=14)
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['HMG','HET'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    plt.xlim(306,326)
    plt.xticks([310,315,320,325],rotation=30)
    plt.ylim(0,5.1)
    
    
    plt.subplot(2,3,i+3)
    q1=flms[day]['qv'][t_i,0:km5l+10]*1000
    qc=flts[day]['qv'][t_i,0:km5l+10]*1000
    plt.plot(q1,alt,':',c='cornflowerblue')
    plt.plot(qc,alt,c='cornflowerblue')
    
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['HMG','HET'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    plt.xlabel('$q_t$ ($g\ kg^{-1}$)',fontsize=14)
    plt.xticks([0,5,10],['','5','10'])
    plt.xlim(-1,13)
    plt.ylim(0,5.1)
    
    i=i+1
plt.subplots_adjust(hspace=.3)
plt.show()    

# %%
7*12

# %%
fcs[day]['frac'][1]

# %% [markdown]
# # Advection Profiles
#

# %%
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####

#'T_f[K\\s]', 'rtm_f[kg\\kg\\s]'

plt.figure(figsize=(6,8))
alt=ft1s[day]['altitude'][0:km5+10]/1000
i=1
t_i=84

for day in core_days:
    plt.subplot(2,3,i)
    t1=frc0[day]['T_f[K\\s]'][t_i,0:km5+10]
    t2=(frcc[day]['T_f[K\\s]'][t_i,0:km5+10]-t1)*fcs[day]['frac'][0]+(frch[day]['T_f[K\\s]'][t_i,0:km5+10]-t1)*fcs[day]['frac'][1]
    plt.plot(t1*10**5,alt,':',c='indianred')
    plt.plot(t2*10**5,alt,'--',c='indianred')
    #plt.plot(tc,alt,c='indianred')
    
    plt.title(labels[i-1],fontsize=14)
    if i ==2:
        plt.xlabel(r'$\frac{d\theta}{dt}$ ($10^{-5}\ K s^{-1}$)',fontsize=14)
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['Large Scale','Circulation'],handlelength=1,loc='upper left')
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    
    plt.xlim(-11,4)
    plt.xticks(rotation=30)
    plt.ylim(0,5.1)
    
    
    plt.subplot(2,3,i+3)
    q1=frc0[day]['rtm_f[kg\\kg\\s]'][t_i,0:km5+10]
    q2=(frcc[day]['rtm_f[kg\\kg\\s]'][t_i,0:km5+10]-q1)*fcs[day]['frac'][0]+(frch[day]['rtm_f[kg\\kg\\s]'][t_i,0:km5+10]-q1)*fcs[day]['frac'][1]
    plt.plot(q1*10**7,alt,':',c='cornflowerblue')
    plt.plot(q2*10**7,alt,'--',c='cornflowerblue')
    #plt.plot(qc,alt,c='cornflowerblue')
    
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['Large Scale','Circulation'],handlelength=1,loc='upper left')
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    if i ==2:
        plt.xlabel(r'$\frac{dq}{dt}$ ($10^{-7}\ g\ kg^{-1} s^{-1}$)',fontsize=14)
    #plt.xticks([0,5,10],['','5','10'])
    plt.xlim(-0.5,4)
    plt.ylim(0,5.1)
    
    i=i+1
    
plt.subplots_adjust(hspace=.4)
plt.show()

# %%

# %% [markdown]
# # 

# %% [markdown]
# # LWP CLUBB ONLY

# %%
#### LWP CLUBB ONLY ####
#### LWP CLUBB ONLY ####
#### LWP CLUBB ONLY ####
fig=plt.figure(figsize=(8,7),dpi=300)
gs = GridSpec(6, 2, figure=fig)
i=0
times=np.linspace(7,23,955)

ytxt=[155,25.5,78]

for day in core_days:
    clr1='dimgrey'
    clr2='grey'
    ccld='cornflowerblue'
    cwrm='indianred'
    
    
    lwp_3=fss[day]['lwp'][:,0,0,0]*1000
    lwp_2=fs2s[day]['lwp'][:,0,0,0]*1000
    lwp_1=fs1s[day]['lwp'][:,0,0,0]*1000
    lwp_3c=fscs[day]['lwp'][:,0,0,0]*1000
    lwp_3h=fshs[day]['lwp'][:,0,0,0]*1000
    lwp_2c=fs2cs[day]['lwp'][:,0,0,0]*1000
    lwp_2h=fs2hs[day]['lwp'][:,0,0,0]*1000
    
    ax = fig.add_subplot(gs[i*2:i*2+2, 0])
    ax.plot(times,lwp_1,':',c=clr1)
    ax.plot(times,lwp_2,'--',c=clr1)
    ax.plot(times,lwp_3,'-',c='black')
    ax.set_xlim(8,20.5)
    ax.set_xticks([8,10,12,14,16,18,20])
    ax.text(8.2,ytxt[i],labels[i],fontsize=15)
    #ax.set_title(labels[i],loc='right')
    #ax.set_ylabel(labels[i],loc='right')
    if i==1:
        ax.set_ylabel('Liquid Water Path (LWP) ($g\ m^{2}$)',fontsize=15)
    if i==0:
        ax.legend(['SC','IC','TCM'],loc='lower left')
    if i==2:
        ax.set_xticklabels(['8:00','','12:00','','16:00','','20:00'])
        ax.set_xlabel('Local Time',fontsize=15)
    else:
        ax.set_xticklabels([])
    
    
    ax=fig.add_subplot(gs[i*2,1])
    ax.plot(times,lwp_2-lwp_1,c=clr2)
    ax.plot(times,lwp_2c-lwp_1,c=ccld,linewidth=1)
    ax.plot(times,lwp_2h-lwp_1,c=cwrm,linewidth=1)
    ax.set_xlim(8,20.5)
    ax.set_xticks([10,12,14,16,18,20])
    ax.set_xticklabels([])
    if i==0:
        ax.legend(['$\Delta$(IC-SC)','Cool','Warm'],handlelength=1,borderpad=.2,labelspacing=.2,loc='upper left')
    
    ax=fig.add_subplot(gs[i*2+1,1])
    ax.plot(times,lwp_3-lwp_1,c=clr2)
    ax.plot(times,lwp_3c-lwp_1,c=ccld,linewidth=1)
    ax.plot(times,lwp_3h-lwp_1,c=cwrm,linewidth=1)
    ax.set_xlim(8,20.5)
    ax.set_xticks([8,10,12,14,16,18,20])
    if i==0:
        ax.legend(['$\Delta$(TCM-SC)','Cool','Warm'],handlelength=1,borderpad=.2,labelspacing=.2,loc='upper left')
    if i==2:
        ax.set_xticklabels(['8:00','','12:00','','16:00','','20:00'])
        ax.set_xlabel('Local Time',fontsize=15)
    else:
        ax.set_xticklabels([])
    i=i+1
plt.show()

# %% [markdown]
# # 

# %% [markdown]
# # VARIANCES LES CLUBB

# %%
#### VARIANCES LES CLUBB ####
#### VARIANCES LES CLUBB ####
#### VARIANCES LES CLUBB ####
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####
#### MEAN PROFILE COMPARISON ####
plt.figure(figsize=(6,8))
alt=flts[day]['z'][50,0:km5l+10]/1000
altc=ft1s[day]['altitude'][0:km5+10]/1000
i=1
tt=16
t_i=6*(tt-7)
t_ic=60*(tt-7)
c1les='olivedrab'
c2les='darkolivegreen'
c1clb='dimgrey'
c2clb='black'
for day in core_days:
    plt.subplot(2,3,i)
    t1=flms[day]['thl2'][t_i,0:km5l+10]
    tc=flts[day]['thl2'][t_i,0:km5l+10]
    
    t1c=fm1s[day]['thlp2'][t_ic,0:km5+10,0,0]
    tcc=fms[day]['thlp2'][t_ic,0:km5+10,0,0]
    
    plt.plot(t1c,altc,':',c=c1clb)
    plt.plot(tcc,altc,c=c2clb)
    plt.plot(t1,alt,':',c=c1les)
    plt.plot(tc,alt,c=c2les)
    
    
    plt.title(labels[i-1],fontsize=14)
    plt.xlabel(r'$\theta_v^2$ ($K^2$)',fontsize=14)
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['SC','TCM','HMG','HET'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    #plt.xlim(310,326)
    #plt.xticks([315,320,325])
    plt.ylim(0,5.1)
    
    
    plt.subplot(2,3,i+3)
    q1=flms[day]['qv2'][t_i,0:km5l+10]*1000
    qc=flts[day]['qv2'][t_i,0:km5l+10]*1000
    
    q1c=fm1s[day]['rtp2'][t_ic,0:km5+10,0,0]*1000
    qcc=fms[day]['rtp2'][t_ic,0:km5+10,0,0]*1000
    
    plt.plot(q1c,altc,':',c=c1clb)
    plt.plot(qcc,altc,c=c2clb)
    plt.plot(q1,alt,':',c=c1les)
    plt.plot(qc,alt,c=c2les)
    
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['SC','TCM','HMG','HET'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    plt.xlabel('$q_t^2$ ($g\ kg^{-1}$)',fontsize=14)
    plt.xticks(rotation=45)
    #plt.xticks([0,5,10],['','5','10'])
    #plt.xlim(-1,13)
    plt.ylim(0,5.1)
    
    i=i+1
plt.subplots_adjust(hspace=.25)
plt.show()    

# %% [markdown]
# # 

# %% [markdown]
# # 2D CLOUD STRUCTURE 

# %%
#### 2D CLOUD STRUCTURE ####
#### 2D CLOUD STRUCTURE ####
#### 2D CLOUD STRUCTURE ####
#'rcm/lwc'
i=0
cmap='Blues'
frc=1.3
vmin=0
vmax=.03
labels=['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)']
sns.set_style("darkgrid", {"axes.facecolor": "#EAEAE2"})

lest=np.linspace(7,22,flts[day]['lwc'][:].shape[0])
clbt=np.linspace(7,23,ft1s[day]['rcm'][:].shape[0])


fig = plt.figure(figsize=(9,7.5))
subfigs = fig.subfigures(3, 3, hspace=0, height_ratios=[1, 1,.55],frameon=False)
# Build Grid
gd=[]
for j,day in enumerate(core_days):
    cldm_c=ft1s[day]['rcm'][:,0:km5,0,0]*1000
    cldt_c=fts[day]['rcm'][:,0:km5,0,0]*1000
    cldm_l=flms[day]['lwc'][:,0:km5l]*1000
    cldt_l=flts[day]['lwc'][:,0:km5l]*1000
    
    lwp_3=fss[day]['lwp'][:,0,0,0]*1000
    lwp_1=fs1s[day]['lwp'][:,0,0,0]*1000
    lwp_m=flms[day]['lwp'][:]*1000
    lwp_t=flts[day]['lwp'][:]*1000
    
    data=np.concatenate((cldm_c.flatten(),cldt_c.flatten(),cldm_l.flatten(),cldt_l.flatten()))
    
    vmax=np.percentile(data,99.5)
    vmin=0
    
    for i in range(2,-1,-1):
        if i<2:
        # grid for top two
            grid=ImageGrid(subfigs[i,j], 111,  # similar to subplot(111)
                 nrows_ncols=(2, 1),
                 axes_pad=0.02,
                 cbar_mode='single',
                 cbar_location='bottom',
                 cbar_pad=.02,
                 cbar_size="10%")
            if i==0:
                data1=cldm_c[0:-60,:]
                data2=cldt_c[0:-60,:]
                cbarlabel=r'CLUBB Liquid Water ($g\ kg^{-1}$)'
                txt=['SC','TCM']
                grid[0].set_title(labels[j],fontsize=15)
            if i==1:
                data1=cldm_l
                data2=cldt_l
                cbarlabel=r'LES Liquid Water ($g\ kg^{-1}$)'
                txt=['HMG','HET']
    
            data1[data1<=.0002]=float('nan')
            data2[data2<=.0002]=float('nan')
            im=grid[0].imshow(data1.T,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
            im=grid[1].imshow(data2.T,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
            
            grid.cbar_axes[0].colorbar(im,label=cbarlabel)
            grid[0].text(8,frc*4+.1,txt[0])
            grid[1].text(8,frc*4+.1,txt[1])
            
            if j==0:
                grid[0].set_yticks([frc,frc*2,frc*3,frc*4],labels=[1,2,3,4])
                grid[1].set_yticks([frc,frc*2,frc*3,frc*4],labels=[1,2,3,4])
                grid[0].set_ylabel('Altitude ($km$)                  ',fontsize=14)
            else:
                grid[0].set_yticks([frc,frc*2,frc*3,frc*4],labels=[])
                grid[1].set_yticks([frc,frc*2,frc*3,frc*4],labels=[])
                
            grid[0].set_xticks([8,12,16,20],labels=[])
            grid[1].set_xticks([8,12,16,20],labels=[])
            
            grid[1].set_zorder((i+5))
            
        else:
            # plotting for the LWP
            axs=subfigs[i,j].subplots(2,1,gridspec_kw={'height_ratios':[.1,1]})
            axs[0].axis('off')
            ax=axs[1]
            ax.plot(clbt,lwp_1,':',c='dimgrey')
            ax.plot(clbt,lwp_3,'-',c='black')
            ax.plot(lest,lwp_m,':',c='olivedrab')
            ax.plot(lest,lwp_t,'-',c='darkolivegreen')
            xticks=np.array([8,12,16,20])
            ax.set_xticks(xticks,labels=['8:00','12:00','16:00','20:00'])
            ax.set_xlim(7,22)
            ax.set_xlabel('Local Time')
            if j==0:
                ax.set_ylabel('LWP ($kgm^{-2}$)',fontsize=14)
                ax.legend(['SC','TCM','HMG','HET'],handlelength=1,borderpad=.2,labelspacing=.2,loc='upper left')

plt.show()

# %% [markdown]
# # 

# %% [markdown]
# # Profile LES Comparison

# %%
plt.clf()
sns.set_theme()
fig=plt.figure(figsize=(12,5),dpi=300)
gs = GridSpec(3, 8, figure=fig,width_ratios=[.8,1,.3,.8,1,.3,.8,1])
i=0
sigma=10
mn=5
wind_scl=np.array([28,28,28])
zlvl={core_days[0]:[mn,mn+66,mn+99],
      core_days[1]:[mn,mn+66,mn+99],
      core_days[2]:[mn,mn+66,mn+99],
     }
pnts={core_days[0]:[[100,260],[350,60]],
      core_days[1]:[[100,260],[350,380]],
      core_days[2]:[[415,415],[150,260]],
     }
w_reso=50
wticks=np.linspace(w_reso,520-w_reso,int((520-2*w_reso)/w_reso)+1)
x,y=np.meshgrid(wticks,wticks)
tm=17+5
symb=['^','.','*','v']
ylabels=['0.17','2.20','3.20']
for day in core_days:
    print(day,flush=True)
    dlhm='/home/tsw35/xTyc_shared/clasp/fr2/fr2_'+day+'_01/'
    dlht='/home/tsw35/xTyc_shared/clasp/fr2/fr2_'+day+'_00/'
    tfile='diag_d01_'+day[0:4]+'-'+day[4:6]+'-'+day[6:8]+'_'+str(tm)+':00:00'
    fpht=nc.Dataset(dlht+tfile,'r')
    fphm=nc.Dataset(dlhm+tfile,'r')
    alt=fpht['AVP_Z'][0,0:km4l]/1000
    uht=fpht['AVV_U'][0,0:km4l,:,:]
    uhm=fphm['AVV_U'][0,0:km4l,:,:]
    vht=fpht['AVV_V'][0,0:km4l,:,:]
    vpt=fpht['AVV_THV'][0,0:km4l,:,:]
    
    cmsk=vpt[5,:,:]<np.percentile(vpt[5,:,:],25)
    hmsk=vpt[5,:,:]>=np.percentile(vpt[5,:,:],75)
    
    prof_c=np.mean(vpt[:,cmsk],axis=1)
    prof_h=np.mean(vpt[:,hmsk],axis=1)
    p1x=pnts[day][0][0]
    p1y=pnts[day][0][1]
    p2x=pnts[day][1][0]
    p2y=pnts[day][1][1]
    #prof_h=np.mean(vpt[:,p1x-5:p1x+5,p1y-5:p1y+5],axis=(1,2))
    #prof_c=np.mean(vpt[:,p2x-5:p2x+5,p2y-5:p2y+5],axis=(1,2))
    
    zs=zlvl[day][::-1]
    
    for j in range(3):
        print(j,flush=True)
        ax1 = fig.add_subplot(gs[j, i*3+1])
        d_u=np.zeros((len(wticks),len(wticks)))
        d_v=np.zeros((len(wticks),len(wticks)))
        z=zs[j]
        vpt_j=sci.ndimage.filters.gaussian_filter(vpt[z,:,:],sigma,mode='reflect')
        datau=sci.ndimage.filters.gaussian_filter(uht[z,:,:]-np.mean(uht[z,:,:]),sigma*2,mode='reflect')
        datav=sci.ndimage.filters.gaussian_filter(vht[z,:,:]-np.mean(vht[z,:,:]),sigma*2,mode='reflect')
        for ii in range(len(wticks)):
            for jj in range(len(wticks)):
                d_u[ii,jj]=datau[int(wticks[ii]),int(wticks[jj])]
                d_v[ii,jj]=datav[int(wticks[ii]),int(wticks[jj])]
        
        im=ax1.imshow(vpt_j,origin='lower',cmap='coolwarm',vmin=np.percentile(vpt_j,2),vmax=np.percentile(vpt_j,98))
        ax1.quiver(x,y,d_u[:,:],d_v[:,:],facecolor='white',edgecolor='black',width=.02,scale=wind_scl[i],alpha=1,headlength=5,linewidth = .1,zorder=2)
        ax1.contour(cmsk,colors='navy',linewidths=0.5,linestyles='dotted',zorder=1)
        ax1.contour(hmsk,colors='maroon',linewidths=0.5,linestyles='dotted',zorder=1)
        #ax1.scatter(50,425,marker=symb[j],s=75,c='black')
        #ax1.scatter(pnts[day][0][0],pnts[day][0][1],marker='o',s=75,c='black')
        #ax1.scatter(pnts[day][1][0],pnts[day][1][1],marker='o',s=75,c='black')
        #ax1.annotate(str(z*30/1000)+' km', (25,425))
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_ylabel(ylabels[2-j]+' km')
        ax1.yaxis.set_label_position("right")
        #plt.colorbar(im,ax=ax1,location='left',label=r'$\theta_v$ (K)')
        
    ax2 = fig.add_subplot(gs[0:4, i*3])
    ax2.plot(prof_c,alt,c='darkblue')
    ax2.plot(prof_h,alt,c='darkred')
    #for j in range(4):
    #    z=zs[j]
    #    ax2.scatter(prof_c[z],alt[z],marker=symb[j],s=75,c='black',zorder=5)
    ax2.set_title(labels[i],fontsize=14)
    #ax2.yaxis.tick_right()
    #ax2.yaxis.set_label_position("right")
    if i==0:
        ax2.set_ylabel('Altitude (km)',fontsize=14)
        ax2.set_yticks(alt[zs],labels=ylabels[::-1],rotation=90)
    else:
        ax2.set_yticks(alt[zs],labels=[])
    ax2.set_ylim(0,4)
    ax2.set_xlabel(r'$\theta_v$ (K)')
    ax2.grid(which='both',axis='y',lw=3,ls='--')
    if i<2:
        ax3 = fig.add_subplot(gs[0:4, i*3+2])
        ax3.axis('off')
    if i==0:
        ax2.legend(['Cool','Warm'],handlelength=1)
    
    i=i+1
plt.subplots_adjust(hspace=0.05,wspace=0)
plt.savefig('junk.png')

# %%

# %%

# %%

# %%
sigma=10
w_reso=40
zlist=[10,40,80,100]
wticks=np.linspace(w_reso,520-w_reso,int((520-2*w_reso)/w_reso)+1)
x,y=np.meshgrid(wticks,wticks)
frames=len(zlist)
d_u = np.zeros((frames,len(wticks),len(wticks)))
d_v = np.zeros((frames,len(wticks),len(wticks)))
d_vpt=np.zeros((frames,520,520))
d_h=np.zeros((frames,))
d_c=np.zeros((frames,))
for k in range(len(zlist)):
    z=zlist[k]
    d_vpt[k,:,:]=sci.ndimage.filters.gaussian_filter(vpt[z,:,:],sigma,mode='reflect')
    d_h[k]=np.mean(d_vpt[k,:,50:150])
    d_c[k]=np.mean(d_vpt[k,:,300:500])
    datau=sci.ndimage.filters.gaussian_filter(uht[z,:,:]-np.mean(uht[z,:,:]),sigma*2,mode='reflect')
    datav=sci.ndimage.filters.gaussian_filter(vht[z,:,:]-np.mean(vht[z,:,:]),sigma*2,mode='reflect')
    for i in range(len(wticks)):
        for j in range(len(wticks)):
            d_u[k,i,j]=datau[int(wticks[i]),int(wticks[j])]
            d_v[k,i,j]=datav[int(wticks[i]),int(wticks[j])]
alt=np.array(zlist)*30

# %%
data=d_vpt[i,:,:]
ax1.clear()
ax1.imshow(data,origin='lower',cmap='coolwarm',vmin=np.mean(data)-.9,vmax=np.mean(data)+.9)
ax1.quiver(x,y,d_u[i,:,:],d_v[i,:,:],color='whitesmoke',width=.010,scale=40,alpha=.95)
ax1.scatter(400,260,s=75,color='steelblue')
ax1.scatter(100,260,s=75,color='firebrick')
ax1.axis('off')
ax1.set_title('Virtual Potential\nTemperature Field')

#plt.colorbar(sc)

#if i==0:
#    fig.colorbar(sc,ax=ax1)

ax2.clear()
ax2.plot(d_c,alt,c='dodgerblue')
ax2.plot(d_h,alt,c='indianred')
ax2.scatter(d_c[i],alt[i],s=75,color='steelblue')
ax2.scatter(d_h[i],alt[i],s=75,color='firebrick')
ax2.set_ylabel('Altitude ($m$)')
ax2.set_xlabel('VPT ($K$)')
ax2.set_xlim(312,320.1)
ax2.set_ylim(-25,3550)
fig.subplots_adjust(wspace=.45, hspace=.45)

# %%

# %%

# %%

# %%
for day in core_days:
    plt.figure(dpi=50)
    data1=np.max(fcs[day]['wbar'][0,:,:],axis=1)
    data2=np.max(fcs[day]['wbar'][1,:,:],axis=1)
    plt.plot(data1)
    plt.plot(data2)
    plt.title(day)
plt.show()

# %%

# %%
#### VARIANCES LES CLUBB ####
#### VARIANCES LES CLUBB ####
#### VARIANCES LES CLUBB ####
plt.figure(figsize=(6,8))
alt=flts[day]['z'][50,0:km5l+10]/1000
altc=ft1s[day]['altitude'][0:km5+10]/1000
i=1
tt=17
t_i=6*(tt-7)
t_ic=60*(tt-7)
c1les='olivedrab'
c2les='darkolivegreen'
c1clb='dimgrey'
c2clb='black'
for day in core_days:
    plt.subplot(2,3,i)
    t1=flms[day]['w2'][t_i,0:km5l+10]#+flms[day]['v2'][t_i,0:km5l+10]
    tc=flts[day]['w2'][t_i,0:km5l+10]#+flts[day]['v2'][t_i,0:km5l+10]
    
    t1c=fm1s[day]['wp2'][t_ic,0:km5+10,0,0]#+fm1s[day]['vp2'][t_ic,0:km5+10,0,0]
    tcc=fms[day]['wp2'][t_ic,0:km5+10,0,0]#+fms[day]['vp2'][t_ic,0:km5+10,0,0]
    
    plt.plot(t1c,altc,':',c=c1clb)
    plt.plot(tcc,altc,c=c2clb)
    plt.plot(t1,alt,':',c=c1les)
    plt.plot(tc,alt,c=c2les)
    
    
    plt.title(labels[i-1],fontsize=14)
    plt.xlabel(r'$\theta_v^2$ ($K^2$)',fontsize=14)
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['SC','TCM','HMG','HET'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    #plt.xlim(310,326)
    #plt.xticks([315,320,325])
    plt.ylim(0,5.1)
    
    
    plt.subplot(2,3,i+3)
    q1=flms[day]['w'][t_i,0:km5l+10]*1000
    qc=flts[day]['w'][t_i,0:km5l+10]*1000
    
    q1c=ft1s[day]['wm'][t_ic,0:km5+10,0,0]*1000
    qcc=fts[day]['wm'][t_ic,0:km5+10,0,0]*1000
    
    plt.plot(q1c,altc,':',c=c1clb)
    plt.plot(qcc,altc,c=c2clb)
    plt.plot(q1,alt,':',c=c1les)
    plt.plot(qc,alt,c=c2les)
    
    if i==1:
        plt.yticks([1,2,3,4,5])
        plt.ylabel('Altitude ($km$)',fontsize=14)
        plt.legend(['SC','TCM','HMG','HET'],handlelength=1)
    else:
        plt.yticks([1,2,3,4,5],labels=[])
    plt.xlabel('$r_t^2$ ($g\ kg^{-1}$)',fontsize=14)
    plt.xticks(rotation=45)
    #plt.xticks([0,5,10],['','5','10'])
    #plt.xlim(-1,13)
    plt.ylim(0,5.1)
    
    i=i+1
plt.subplots_adjust(hspace=.25)
plt.show()    

# %%
i=1
for day in core_day:
    plt.subplot(3,2,1)
    plt.imshow(data1.T,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)

# %%
#### LWP CLUBB ONLY ####
#### LWP CLUBB ONLY ####
#### LWP CLUBB ONLY ####
fig=plt.figure(figsize=(8,7),dpi=300)
gs = GridSpec(6, 2, figure=fig)
i=0
times=np.linspace(7,23,955)

ytxt=[155,25.5,78]

for day in core_days:
    clr1='dimgrey'
    clr2='grey'
    ccld='cornflowerblue'
    cwrm='indianred'
    
    
    lwp_3=fss[day]['lwp'][:,0,0,0]*1000
    lwp_2=fs2s[day]['lwp'][:,0,0,0]*1000
    lwp_1=fs1s[day]['lwp'][:,0,0,0]*1000
    lwp_3c=fscs[day]['lwp'][:,0,0,0]*1000
    lwp_3h=fshs[day]['lwp'][:,0,0,0]*1000
    lwp_2c=fs2cs[day]['lwp'][:,0,0,0]*1000
    lwp_2h=fs2hs[day]['lwp'][:,0,0,0]*1000
    
    ax = fig.add_subplot(gs[i*2:i*2+2, 0])
    ax.plot(times,lwp_1,':',c=clr1)
    ax.plot(times,lwp_2,'--',c=clr1)
    ax.plot(times,lwp_3,'-',c='black')
    ax.set_xlim(8,20.5)
    ax.set_xticks([8,10,12,14,16,18,20])
    #ax.text(8.2,ytxt[i],labels[i],fontsize=15)
    #ax.set_title(labels[i],loc='right')
    #ax.set_ylabel(labels[i],loc='right')
    if i==1:
        ax.set_ylabel('Liquid Water Path (LWP) ($g\ m^{2}$)',fontsize=15)
    if i==0:
        ax.legend(['SC','IC','TCM'],loc='upper left')
    if i==2:
        ax.set_xticklabels(['8:00','','12:00','','16:00','','20:00'])
        ax.set_xlabel('Local Time',fontsize=15)
    else:
        ax.set_xticklabels([])
    
    ax=fig.add_subplot(gs[i*2:i*2+2,1])
    ax.plot(times,lwp_3-lwp_1,c=clr2)
    ax.plot(times,lwp_3c-lwp_1,c=ccld,linewidth=1)
    ax.plot(times,lwp_3h-lwp_1,c=cwrm,linewidth=1)
    ax.set_xlim(8,20.5)
    ax.set_xticks([8,10,12,14,16,18,20])
    if i==0:
        ax.legend(['$\Delta$(TCM-SC)','Warm','Cool'],handlelength=1,borderpad=.2,labelspacing=.2,loc='upper left')
    if i==2:
        ax.set_xticklabels(['8:00','','12:00','','16:00','','20:00'])
        ax.set_xlabel('Local Time',fontsize=15)
    else:
        ax.set_xticklabels([])
    i=i+1
plt.show()

# %%
#### CIRCULATION COMPARISON LES and CLUBB ####
#### CIRCULATION COMPARISON LES and CLUBB ####
#### CIRCULATION COMPARISON LES and CLUBB ####
km5=83
km5l=166
frc=1.25
cc1=.85
vmax=4
vmin=0
testdir=clubb_dir+subdir+'_cpl/'
cmap='terrain'
i=1

plt.figure(figsize=(10,4),dpi=300)

#filelist=['20160625.nc','20170716.nc','20170717.nc','20190707.nc','20180709.nc','20180707.nc','20170802.nc','20150801.nc','20160610.nc','20160716.nc']
for day in ['20160625','20170716','20170717']:
    file='sgp_'+day
        
    day=file[4:]
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    ft=nc.Dataset(testdir+file+'/k_2/agg_outzt.nc','r')
    
    fp=nc.Dataset(les_prmd+day+'.nc','r')

    runpm=open(testdir+file+'/tw_run_param.txt','r').readlines()
    
    for ln in runpm:
        if 'c_ur' in ln:
            c_ur=re.findall('\d+\.\d+',ln)[0]
            
    ur0 =fc['u_r0'][:,:km5]
    uu  =ft['um'][:,:km5,0,0]
    vv  =ft['vm'][:,:km5,0,0]
    wfrac=fc['W_frac'][:]
    ur=fc['u_r'][:,:km5]
    ur[ur==-1]=0#float('nan')
    ur[ur==0]=0#float('nan')
    ur=ur[0:-12,:]
    maxzh=fc['z_circh'][:-12]
    data= np.abs(ur.T)
    if day=='20160625':
        data=data/fc['W_frac'][0]
    
    data=circ2thick(data,maxzh/60,day=='20170717')
    
    plt.subplot(2,3,i)
    #data[data==0]=float('nan')
    plt.imshow(data,cmap=cmap,origin='lower',vmin=vmin,vmax=vmax,extent=(7,22,0,5*frc))
    #if i==3:
    #    plt.colorbar(label=r'$u_r$ ($ms^{-1}$)')
    #plt.title(dayst[day])
    #plt.title(dayst[day]+'\n'+labels[i-1]+'                               ',fontsize=15)
    
    #plt.plot(np.linspace(7,22,len(maxzh)),maxzh/1000*frc,'k--')
    plt.ylim(0,frc*5)
    plt.xlim(8,20)
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([9,12,15,18],labels=[])
    if i==1:
        plt.ylabel('Altitude ($km$)                        ',fontsize=15)
    
    #plt.xlabel('Local Time')
    
    
    
    
    
    plt.subplot(2,3,i+3)
    u90=fps[day]['u90'][:]
    u10=fps[day]['u10'][:]
    v90=fps[day]['v90'][:]
    v10=fps[day]['v10'][:]
    
    vpt=fps[day]['vpt'][:]
    
    #plt.subplot(1,2,1)
    #plt.title(dayst[day])
    #data=u90[:,0:166].T-u10[:,0:166].T
    data=np.sqrt(((u90[:,0:166].T-u10[:,0:166].T)/2)**2+((v90[:,0:166].T-v10[:,0:166].T)/2)**2)
    
    nch1=(u90[:,0:km5l].T/np.abs(u90[:,0:km5l].T))==(u10[:,0:km5l].T/np.abs(u10[:,0:km5l].T))
    nch2=(v90[:,0:km5l].T/np.abs(v90[:,0:km5l].T))==(v10[:,0:km5l].T/np.abs(v10[:,0:km5l].T))
    
    
    '''
    data[nch1&nch2]=float('nan')
    
    if day=='20160625':
        data[int(3.5/5*km5l):,:]=float('nan')
        
    if day=='20170716':
        data[int(3/5*km5l):,:]=float('nan')
        
    if day=='20170717':
        data[int(3.1/5*km5l):,:]=float('nan')
    '''
    
    #data=data/np.max(data,axis=0)
    #data[data<1]=float('nan')
    plt.imshow(data,cmap=cmap,origin='lower',vmin=vmin,vmax=vmax,extent=(7,22,0,5*frc))
    #plt.colorbar(label=r'$u_r$ ($ms^{-1}$)')
    #if i==3:
    #    plt.colorbar(label=r'$u_r$ ($ms^{-1}$)')
    
    #plt.plot(np.linspace(7,22,16),maxzh*30*frc/1000,'k--')
    #plt.plot(np.linspace(7,22,16),maxzl*30*frc/1000,'k-')
    #plt.plot(np.linspace(7,22,16),minzh*30*frc/1000,'k:')
    plt.ylim(0,frc*5)
    plt.xlim(8,20)
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([9,12,15,18],labels=['9:00','12:00','15:00','18:00'])
    if i==2:
        plt.xlabel('Local Time',fontsize=15)
    #plt.legend([r'$z_{max_1}$'],loc='upper left',framealpha=.9)
    
    i=i+1
plt.show()    

# %%
for day in core_days:
    fig=plt.figure(figsize=(4,3))
    ax=fig.add_subplot(111)
    lwp_3=fss[day]['lwp'][:,0,0,0]*1000
    lwp_1=fs1s[day]['lwp'][:,0,0,0]*1000
    lwp_m=flms[day]['lwp'][:]*1000
    lwp_t=flts[day]['lwp'][:]*1000
    ax.plot(clbt,lwp_1,':',c='dimgrey')
    ax.plot(clbt,lwp_3,'-',c='black')
    ax.plot(lest,lwp_m,':',c='olivedrab')
    ax.plot(lest,lwp_t,'-',c='darkolivegreen')
    xticks=np.array([8,12,16,20])
    ax.set_xticks(xticks,labels=['8:00','12:00','16:00','20:00'])
    ax.set_xlim(7,22)
    ax.set_xlabel('Local Time')
    ax.set_ylabel('LWP ($kgm^{-2}$)',fontsize=14)
    ax.legend(['SC','TCM','HMG','HET'],handlelength=1,borderpad=.2,labelspacing=.2,loc='upper left')
plt.show()

# %%
print()

# %%
alt=fths[day]['altitude'][:]
ttt=7*60
plt.plot(ftcs[day]['thvm'][ttt,:,0,0],alt)
plt.plot(fths[day]['thvm'][ttt,:,0,0],alt)
plt.show()

# %%
data=np.abs(fts[day]['vm'][:,0:km5,0,0].T)
plt.imshow(data,cmap='terrain',origin='lower',extent=(7,22,0,5*frc))
plt.colorbar()
plt.show()

# %%
fcs[day].variables

# %%
fp=nc.Dataset(clubb_dir+'test_cpl/k_2/')


# %%
def index_time(date,clb_st=7,les_st=7,var_st=datetime.datetime(2012,5,1,0,0)):
    clb_i=(date.hour-clb_st)*60+date.minute
    if clb_i<0:
        clb_i=clb_i+60*24
    les_i=int((date.hour-les_st)*6+date.minute/10)
    var_i=int((date-var_st).total_seconds()/60/60)+clb_st
    return clb_i,les_i,var_i


# %%
varanal='/home/tsw35/soteria/clubb/data/sgp60varanarap_2012-2019.nc'
fpv=nc.Dataset(varanal,'r')
var_p=fpv['lev'][:]

# %%
x,y,var_i1=index_time(datetime.datetime(2016,6,25,0,0))
x,y,var_i2=index_time(datetime.datetime(2017,7,17,0,0))
x,y,var_i3=index_time(datetime.datetime(2018,7,9,0,0))
var_i=[var_i1,var_i2,var_i3]

# %%
#### 2D CLOUD STRUCTURE ####
#### 2D CLOUD STRUCTURE ####
#### 2D CLOUD STRUCTURE ####
#'rcm/lwc'
i=0
cmap='Blues'
frc=1.3
vmin=0
vmax=.03
labels=['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)']
sns.set_style("darkgrid", {"axes.facecolor": "#EAEAE2"})

lest=np.linspace(7,22,flts[day]['lwc'][:].shape[0])
clbt=np.linspace(7,23,ft1s[day]['rcm'][:].shape[0])
vart=np.linspace(7,22,16)


fig = plt.figure(figsize=(9,7.5))
subfigs = fig.subfigures(3, 3, hspace=0, height_ratios=[1, 1,.55],frameon=False)
# Build Grid
gd=[]

for j,day in enumerate(core_days):
    cldm_c=ft1s[day]['rcm'][:,0:km5,0,0]*1000
    cldt_c=fts[day]['rcm'][:,0:km5,0,0]*1000
    cldm_l=flms[day]['lwc'][:,0:km5l]*1000
    cldt_l=flts[day]['lwc'][:,0:km5l]*1000
    
    lwp_3=fss[day]['lwp'][:,0,0,0]*1000
    lwp_1=fs1s[day]['lwp'][:,0,0,0]*1000
    lwp_m=flms[day]['lwp'][:]*1000
    lwp_t=flts[day]['lwp'][:]*1000
    lwp_v=fpv['LWP'][var_i[j]:var_i[j]+16]*10000
    
    data=np.concatenate((cldm_c.flatten(),cldt_c.flatten(),cldm_l.flatten(),cldt_l.flatten()))
    
    vmax=np.percentile(data,99.5)
    vmin=0
    
    for i in range(2,-1,-1):
        if i<2:
        # grid for top two
            grid=ImageGrid(subfigs[i,j], 111,  # similar to subplot(111)
                 nrows_ncols=(2, 1),
                 axes_pad=0.02,
                 cbar_mode='single',
                 cbar_location='bottom',
                 cbar_pad=.02,
                 cbar_size="10%")
            if i==0:
                data1=cldm_c[0:-60,:]
                data2=cldt_c[0:-60,:]
                cbarlabel=r'CLUBB Liquid Water ($g\ kg^{-1}$)'
                txt=['SC','TCM']
                grid[0].set_title(labels[j],fontsize=15)
            if i==1:
                data1=cldm_l
                data2=cldt_l
                cbarlabel=r'LES Liquid Water ($g\ kg^{-1}$)'
                txt=['HMG','HET']
    
            data1[data1<=.0002]=float('nan')
            data2[data2<=.0002]=float('nan')
            im=grid[0].imshow(data1.T,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
            im=grid[1].imshow(data2.T,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
            
            grid.cbar_axes[0].colorbar(im,label=cbarlabel)
            grid[0].text(8,frc*4+.1,txt[0])
            grid[1].text(8,frc*4+.1,txt[1])
            
            if j==0:
                grid[0].set_yticks([frc,frc*2,frc*3,frc*4],labels=[1,2,3,4])
                grid[1].set_yticks([frc,frc*2,frc*3,frc*4],labels=[1,2,3,4])
                grid[0].set_ylabel('Altitude ($km$)                  ',fontsize=14)
            else:
                grid[0].set_yticks([frc,frc*2,frc*3,frc*4],labels=[])
                grid[1].set_yticks([frc,frc*2,frc*3,frc*4],labels=[])
                
            grid[0].set_xticks([8,12,16,20],labels=[])
            grid[1].set_xticks([8,12,16,20],labels=[])
            
            grid[1].set_zorder((i+5))
            
        else:
            # plotting for the LWP
            axs=subfigs[i,j].subplots(2,1,gridspec_kw={'height_ratios':[.1,1]})
            axs[0].axis('off')
            ax=axs[1]
            ax.plot(clbt,lwp_1,':',c='dimgrey')
            ax.plot(clbt,lwp_3,'-',c='black')
            ax.plot(lest,lwp_m,':',c='olivedrab')
            ax.plot(lest,lwp_t,'-',c='darkolivegreen')
            ax.plot(vart,lwp_v,'-')
            xticks=np.array([8,12,16,20])
            ax.set_xticks(xticks,labels=['8:00','12:00','16:00','20:00'])
            ax.set_xlim(7,22)
            ax.set_xlabel('Local Time')
            if j==0:
                ax.set_ylabel('LWP ($kgm^{-2}$)',fontsize=14)
                ax.legend(['SC','TCM','HMG','HET'],handlelength=1,borderpad=.2,labelspacing=.2,loc='upper left')

plt.show()

# %%
plt.plot(vart,lwp_v)

# %%
for v in fpv.variables:
    print(v)

# %%
fpv['LWP']

# %%
print(datetime.datetime(2012,5,1,0,0)+datetime.timedelta(seconds=int(fpv['time'][var_i[0]+1])))

# %%
lest

# %%
km4l

# %%
##################################
#### HORIZONTAL CROSS-SECTION ####
##################################
plt.figure(figsize=(9,6))
i=0
dlht='/home/tsw35/xTyc_shared/clasp/fr2/'
crx=[400,400,400]

for j,day in enumerate(core_days):
    fls3d=nc.Dataset(dlht+'fr2_'+day+'_00/'+'diag_d01_'+day[0:4]+'-'+day[4:6]+'-'+day[6:8]+'_20:00:00','r')
    plt.subplot(3,1,j+1)
    data_h=np.zeros((520,km4l+20))
    data_u=np.zeros((520,km4l+20))
    for i in range(20):
        data_h[:,i]=fls3d['AVV_THV'][0,10,crx[j],:]
        data_u[:,i]=float('nan')
    for i in range(km4l):
        data_u[:,i+20]=fls3d['AVV_U'][0,i,crx[j],:]
        data_h[:,i+20]=float('nan')
    data_u=data_u-np.nanmean(data_u[:,20:53])
    vmx=np.nanpercentile(np.abs(data_u),97)
    plt.imshow(data_h.T,origin='lower',cmap='coolwarm')
    plt.imshow(data_u.T,origin='lower',cmap='PRGn',vmin=-5,vmax=5)
    plt.ylabel('   Altitude (km)',fontsize=14)
    plt.colorbar(label='West-East Velocity')
    plt.xticks([60,140,220,300,380,460],[0,20,40,60,80,100])
    yticks=[20,20+33,20+66,20+100,20+km4l,]
    yticklabel=[0,1,2,3,4]
    plt.yticks(yticks,yticklabel)
    if j==2:
        plt.xlabel('West-East Distance (km)')

# %%
##################################
#### HORIZONTAL CROSS-SECTION ####
##################################
plt.figure(figsize=(9,6))
i=0
dlht='/home/tsw35/xTyc_shared/clasp/fr2/'
crx=[400,400,400]

for j,day in enumerate(core_days):
    fls3d=nc.Dataset(dlht+'fr2_'+day+'_00/'+'diag_d01_'+day[0:4]+'-'+day[4:6]+'-'+day[6:8]+'_20:00:00','r')
    plt.subplot(3,1,j+1)
    data_h=np.zeros((520,km4l+20))
    data_u=np.zeros((520,km4l+20))
    for i in range(20):
        data_h[:,i]=fls3d['AVV_THV'][0,10,crx[j],:]
        data_u[:,i]=float('nan')
    for i in range(km4l):
        data_u[:,i+20]=fls3d['AVV_U'][0,i,crx[j],:]
        data_h[:,i+20]=float('nan')
    data_u=data_u-np.nanmean(data_u[:,20:53])
    data_u[:,20:]=sci.ndimage.filters.gaussian_filter(data_u[:,20:],5,mode='reflect')
    vmx=np.nanpercentile(np.abs(data_u),97)
    plt.imshow(data_h.T,origin='lower',cmap='coolwarm')
    plt.imshow(data_u.T,origin='lower',cmap='PRGn',vmin=-5,vmax=5)
    plt.ylabel('   Altitude (km)',fontsize=14)
    plt.colorbar(label='West-East Velocity')
    plt.xticks([60,140,220,300,380,460],[0,20,40,60,80,100])
    yticks=[20,20+33,20+66,20+100,20+km4l,]
    yticklabel=[0,1,2,3,4]
    plt.yticks(yticks,yticklabel)
    if j==2:
        plt.xlabel('West-East Distance (km)')

# %%
sci.ndimage.filters.gaussian_filter(uht[z,:,:]-np.mean(uht[z,:,:]),sigma*2,mode='reflect')

# %%
fls3d['AVV_U']

# %%
for j,day in enumerate(core_days):
    fls3d=nc.Dataset(dlht+'fr2_'+day+'_00/'+'diag_d01_'+day[0:4]+'-'+day[4:6]+'-'+day[6:8]+'_20:00:00','r')
    plt.subplot(3,1,j+1)
    plt.imshow(fls3d['AVV_THV'][0,10,200:,:],cmap='coolwarm')
    

# %%

# %%
7*60

# %%

for day in core_days:
    dataTc=ftcs[day]['T_in_K'][420,:,0,:]
    dataTVc=ftcs[day]['thvm'][420,:,0,:]
    dataTh=fths[day]['T_in_K'][420,:,0,:]
    dataTVh=ftcs[day]['thvm'][420,:,0,:]
    plt.plot(dataTh)

# %%
day='20170717'

# %%
ft['T_in_K']

# %%
dataTc=ftcs[day]['T_in_K'][420,:,0,:]
dataTVc=ftcs[day]['thvm'][420,:,0,:]
dataTh=fths[day]['T_in_K'][420,:,0,:]
dataTVh=fths[day]['thvm'][420,:,0,:]
tsfh=fshs[day]['T_sfc'][420,0,0,0]
tsfc=fscs[day]['T_sfc'][420,0,0,0]


# %%
plt.plot(dataTVh[0:93],altc)
#plt.plot([0,500],[2.75,2.75],alpha=.25)
plt.plot([dataTVh[0,0],dataTVh[0,0]],[-1,4],alpha=.25)
plt.xlim(310,318)
plt.ylim(0,3.1)

# %%
dataTh[26,0]

# %%
np.where(altc>1.75)[0][0]

# %%
np.where(altc>=3)[0][0]

# %%
(dataTh[51,0]-dataTh[30,0])

# %%
dataTVh[]

# %%
(dataTVh[51,0]-dataTVh[0,0])

# %%
(dataTVh[0,0]-dataTVc[0,0])*15/3

# %%
hclst=fcs[day]['cluster'][7,:,:]==1
cclst=fcs[day]['cluster'][7,:,:]==0

# %%
tsfh=np.mean(fcs[day]['tskin'][7,:][hclst])
tsfh=np.mean(fcs[day]['tskin'][7,:][cclst])

# %%
np.std(fcs[day]['tskin'][7])*15.15

# %%

# %%
fcs[day]['tskin'][7,hclst]

# %%
tsfh

# %%
(dataTVh[51,0]-dataTVh[30,0])*(altc[51]-altc[0])


# %%
def adv_L(wd_,Axx_):
    x1=Axx_[0]/wd_
    x2=Axx_[1]/wd_
    return (x1+x2)/2


# %%
20,5000
wd=fcs[day]['W'][0,0]
wf=fcs[day]['W_frac'][:]
L= adv_L(wd,wf*20*20*5000*5000)
print(L)

# %%
fcs[day].variables

# %%
altc[1]-altc[0]


# %%
def circ_h2(vpt_,L,cc1=1):
    maxzh=np.zeros((16,))
    minzh=np.zeros((16,))
    maxzl=np.zeros((16,))
    print(L)
    for t in range(1,16):
        vpthg=np.gradient(vpt_[t,1,:])
        vpthg[0:5]=0
        cdT=(vpt_[t,1,0]-vpt_[t,0,0])*L/altl
        
        
        vpthgc=np.cumsum(vpthg)
        maxzh[t]=np.argmin(np.abs(vpthgc-cdT))
        print(cdT[int(maxzh[t])])
        print(altl[int(maxzh[t])])
        print()
        vpt_hm=vpt[t,1,int(maxzh[t])]
        maxzl[t]=np.argmax((vpt_[t,0,:]-vpt_hm)>0)-1
        try:
            minzh[t]=np.where((vpt_[t,1,:]-vpt_[t,0,:])<0)[0][0]-1
        except:
            minzh[t]=0
    print()
    return maxzh,maxzl,minzh

# %%
cc1=.85
frc=1.5
km4l=133
fig=plt.figure(figsize=(12,7),dpi=300)
gs = GridSpec(6, 5, figure=fig)
i=0
labels=['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)']

for day in core_days:
    file=day+'.nc'
    u90=fps[day]['u90'][:]
    u10=fps[day]['u10'][:]
    v90=fps[day]['v90'][:]
    v10=fps[day]['v10'][:]
    wd=fcs[day]['W'][0,0]
    wf=fcs[day]['W_frac'][:]
    L= adv_L(wd,wf*20*20*5000*5000)
    
    altc=ft1s[day]['altitude'][:]
    
    vpt=fps[day]['vpt'][:]
    dtsk=2*tsk[day].T
    
    ax1 = fig.add_subplot(gs[i*2:i*2+2, 0:3])
    
    maxzh,maxzl,minzh = circ_h2(vpt,L)

    data=np.sqrt(((u90[:,0:km4l].T-u10[:,0:km4l].T)/2)**2+((v90[:,0:km4l].T-v10[:,0:km4l].T)/2)**2)
    
    data=data/np.max(data,axis=0)
    
    im=ax1.imshow(data,cmap='terrain',origin='lower',extent=(7,22,0,frc*4))
    fig.colorbar(im,ax=ax1,label=r'Normalized $u_r$')
    
    ax1.plot(np.linspace(7,22,16),maxzh*30*frc/1000,'k--')
    ax1.plot(np.linspace(7,22,16),maxzl*30*frc/1000,'k-')
    ax1.plot(np.linspace(7,22,16),minzh*30*frc/1000,'k:')
    ax1.set_ylim(0,frc*4)
    ax1.set_xlim(8,20)
    ax1.set_title(labels[i],loc='left')
    ax1.set_yticks([frc,frc*2,frc*3,frc*4])
    ax1.set_yticklabels([1,2,3,4])
    ax1.set_xticks([8,12,16,20])
    if i==2:
        ax1.set_xticklabels(['8:00','12:00','16:00','20:00'])
        ax1.set_xlabel('Local Time')
    else:
        ax1.set_xticklabels([])
    if i == 1:
        ax1.set_ylabel('Altitude ($km$)',fontsize=16)
    if i == 0:
        ax1.legend([r'$z_{max_1}$',r'$z_{max_2}$',r'$z_{crit}$'],handlelength=1,handletextpad=.6,loc='upper left',framealpha=.9)
    
    i=i+1
ax=fig.add_subplot(gs[1:-1,3:5])
ax.plot(np.linspace(0,7),np.linspace(0,7),'--',color='grey',zorder=0)
ax.scatter(ums[msk]*1.36,u2s[msk],s=6,alpha=.1,zorder=1)
ax.set_xlim(0,7)
ax.set_title(labels[i],loc='left')
ax.set_ylim(0,7)
ax.set_xlabel('Model Velocity ($ms^{-1}$)',fontsize=16)
ax.set_ylabel('LES Velocity ($ms^{-1}$)',fontsize=16)
fig.subplots_adjust(wspace=.35,hspace=.4)
plt.show()

# %%
maxzh

# %%
fps[day].variables

# %%

# %%
for v in fv.variables:
    print(v)

# %%
