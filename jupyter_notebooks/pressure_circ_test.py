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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import datetime
import os
import seaborn as sns
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
sns.set_theme()
plt.rcParams.update({'figure.max_open_warning': 0})

# %%
clubb_dir='/home/tsw35/tyche/clubb/'

# %%

# %%
days=    [20160625,20160716,20160719,20160720,20170609,
          20170626,20170627,20170629,20170705,20170709,
          20170712,20170716,20170717,20170719,20170720,
          20170728,20170826,20170922,20170923,20170924,
          20180522,20180530,20180618,20180619,20180704,
          20180705,20180523,20180707,20180709,20180710,
          20180711,20180712,20180809,20180811,20180916,
          20180917,20190707,20190709,20190714]

# %%
days_=os.listdir(clubb_dir+'sgp_2c_d/')
days_.sort()
daysa=[]
for day in days_:
    if '20190804' in day:
        break
    daysa.append(day[4:])

# %%
urs=[]
for day in daysa:
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d4/sgp_'+str(day)+'/'
    fpc=nc.Dataset(dircp+'k_2/clusters.nc','r')
    urss=fpc['u_r'][:]
    urss[urss<0]=float('nan')
    urs.append(np.nanmax(urss))
#plt.hist(urs,bins=[0,.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7],histtype="step")
plt.hist(urs,bins=[0,.25,.5,.75,1,1.25,1.5,1.75,2,2.5,3,3.5],histtype="step")
#plt.hist(urs,bins=[0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5],histtype="step")
plt.title('')

# %%
urs=[]
for day in daysa:
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d/sgp_'+str(day)+'/'
    fs=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
    break

# %%
plt.plot(fpc['tskin'][:,0,0])

# %%
fpc

# %%
i=0
for day in daysa:
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d/sgp_'+str(day)+'/'
    frc_cc= read_forcings(dircp+'k_2/c_1/arm_forcings_f.in')
    frc_ch= read_forcings(dircp+'k_2/c_2/arm_forcings_f.in')
    frc_0 = read_forcings(dir2c+'k_2/c_1/arm_forcings_f.in')
    
    fscp_c=nc.Dataset(dircp+'k_2/c_1/output/arm_sfc.nc','r')
    fscp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_sfc.nc','r')
    fscp=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
    fs2c_c=nc.Dataset(dir2c+'k_2/c_1/output/arm_sfc.nc','r')
    fs2c_h=nc.Dataset(dir2c+'k_2/c_2/output/arm_sfc.nc','r')
    fs2c=nc.Dataset(dir2c+'k_2/agg_outsfc.nc','r')
    fs1c=nc.Dataset(dir1c+'k_1/c_1/output/arm_sfc.nc','r')
    ftcp_c=nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
    ftcp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
    
    fig=plt.figure(figsize=(10,5))
    
    ax=plt.subplot(2,2,1)
    tfrc=frc_ch['T_f[K\s]'][24:-12,0:150]+frc_cc['T_f[K\s]'][24:-12:,0:150]-2*frc_0['T_f[K\s]'][24:-12:,0:150]
    abmax=.0001#np.percentile(np.abs(tfrc),99)
    im=plt.imshow(tfrc.T,origin='lower',cmap='coolwarm',vmin=-abmax,vmax=abmax,extent=(7,21,0,6))
    #plt.title('T_forcing -- '+str(day))
    plt.ylabel('Elevation (km)')
    plt.xlabel('Time (hr)')
    cax = fig.add_axes([ax.get_position().x1-0.04,ax.get_position().y0,0.02,ax.get_position().height])
    cbar=plt.colorbar(im,cax=cax,label='$\Delta$T K/s')
    cbar.formatter.set_powerlimits((0, 0))
    
    ax=plt.subplot(2,2,2)
    tfrc=frc_ch['rtm_f[kg\kg\s]'][24:-12,0:150]+frc_cc['rtm_f[kg\kg\s]'][24:-12:,0:150]-2*frc_0['rtm_f[kg\kg\s]'][24:-12:,0:150]
    abmax=1.5*10**(-7)#np.percentile(np.abs(tfrc),99)
    im=plt.imshow(tfrc.T,origin='lower',cmap='coolwarm',vmin=-abmax,vmax=abmax,extent=(7,21,0,6))
    #plt.title('RTM_forcing -- '+str(day))
    plt.xlabel('Time (hr)')
    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    cbar=plt.colorbar(im,cax=cax,label='$\Delta$rtm kg/kg/s')
    cbar.formatter.set_powerlimits((0, 0))
    #plt.xlabel('Hour')
    
    '''
    ax=plt.subplot(2,2,3)
    data=ftcp_h['p_in_Pa'][120:-120,0:150,0,0]-ftcp_c['p_in_Pa'][120:-120,0:150,0,0]
    abmax=np.percentile(np.abs(data),95)
    im=plt.imshow(data.T,origin='lower',cmap='coolwarm',vmin=-abmax,vmax=abmax,extent=(7,20,0,6))
    plt.title('Pressure Differences -- '+str(day))
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im,cax=cax)
    
    vpt_h=ftcp_h['thvm'][120:-120,0:150,0,0]
    #vpt_c=ftcp_c['thvm'][120:-120,0:150,0,0]
    maxc=np.zeros((vpt_h.shape[0],))
    for t in range(vpt_h.shape[0]):
        vptg =np.gradient(vpt_h[t,:],)
        vptg[0:5]=0
        vptg=np.cumsum(vptg)
        maxc[t]=np.argmin(np.abs(vptg-.5))/150*6
    
    time=np.linspace(7,20,vpt_h.shape[0])
    ax.plot(time,maxc,c='black')
    '''
    
    lwp_1c=fs1c['lwp'][120:-60,0,0,0]
    time=np.linspace(7,21,lwp_1c.shape[0])
    plt.subplot(2,2,3)
    plt.plot(time,lwp_1c)
    plt.plot(time,fs2c['lwp'][120:-60,0,0,0])
    plt.plot(time,fscp['lwp'][120:-60,0,0,0])
    plt.xlim(7,21)
    plt.xlabel('Time (hr)')
    plt.legend(['1C','2C','CPL'])
    #plt.title('LWP')
    plt.ylabel('LWP')
    
    plt.subplot(2,2,4)
    plt.plot(time,fs2c_c['lwp'][120:-60,0,0,0]-lwp_1c,'b--')
    plt.plot(time,fs2c_h['lwp'][120:-60,0,0,0]-lwp_1c,'r--')
    plt.plot(time,fscp_c['lwp'][120:-60,0,0,0]-lwp_1c,'b-')
    plt.plot(time,fscp_h['lwp'][120:-60,0,0,0]-lwp_1c,'r-')
    plt.xlim(7,21)
    plt.xlabel('Time (hr)')
    plt.legend(['2C wet','2C dry','CPL wet','CPL dry'])
    #plt.title('$\Delta$ LWP from 1C')
    plt.ylabel('$\Delta$ LWP')
    
    i = i+1
    plt.subplots_adjust(hspace=.35,wspace=.45)

# %%
for day in days:
    dir2c = clubb_dir+'sgp_2c_p/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_p/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_p/sgp_'+str(day)+'/'

    ftcp_c=nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
    ftcp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
    fscp_c=nc.Dataset(dircp+'k_2/c_1/output/arm_sfc.nc','r')
    fscp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_sfc.nc','r')
    fmcp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_zm.nc','r')
    fscp=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
    
    Hdif=np.mean(fscp_h['sh'][:,0,0,0]-fscp_c['sh'][:,0,0,0])
    Hmen=np.mean(fscp_h['sh'][:,0,0,0]+fscp_c['sh'][:,0,0,0])/2
    
    fig=plt.figure(figsize=(10,5))

    data=ftcp_h['p_in_Pa'][120:-120,0:150,0,0]-ftcp_c['p_in_Pa'][120:-120,0:150,0,0]
    abmax=10#np.percentile(np.abs(data),95)
    im=plt.imshow(data.T,origin='lower',cmap='coolwarm',vmin=-abmax,vmax=abmax,extent=(7,20,0,6))
    plt.title('Pressure Differences -- '+str(day)+' '+str(int(Hdif))+' '+str(int(Hmen)))
    plt.colorbar()
    
    vpt_h=ftcp_h['thvm'][120:-120,0:150,0,0]
    rig=fmcp_h['Richardson_num'][120:-120,0:150,0,0]
    # vpt gradient cumulative cuttoff 
    plt.plot(time,calc_maxc(vpt_h,1),'k',lineWidth=1)

    # vpt gradient max
    plt.plot(time,calc_maxg(vpt_h,.5),'k--',lineWidth=1)
    
    # richardson number height
    plt.plot(time,calc_blh(rig,.75),'k-.',lineWidth=1)

    
    
    #plt.plot(time,calc_maxc(vpt_h,1),'k',lineWidth=1)
    #plt.plot(time,calc_maxc(vpt_h,1.5),'k',lineWidth=1)
    
    time=np.linspace(7,20,vpt_h.shape[0])
    

# %%
from matplotlib import animation
from IPython.display import HTML
day=20180707
dir2c = clubb_dir+'sgp_2c_p/sgp_'+str(day)+'/'
dir1c = clubb_dir+'sgp_1c_p/sgp_'+str(day)+'/'
dircp = clubb_dir+'sgp_cpl_p/sgp_'+str(day)+'/'

ftcp_c=nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
ftcp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
ft2c_c=nc.Dataset(dir2c+'k_2/c_1/output/arm_zt.nc','r')
ft2c_h=nc.Dataset(dir2c+'k_2/c_2/output/arm_zt.nc','r')

rtm_lo=ft2c_c['rtm'][120:-120,0:150,0,0][times.astype(int),:]
rtm_hi=ftcp_c['rtm'][120:-120,0:150,0,0][times.astype(int),:]
T_lo=ft2c_c['thvm'][120:-120,0:150,0,0][times.astype(int),:]
T_hi=ftcp_c['thvm'][120:-120,0:150,0,0][times.astype(int),:]

from matplotlib import animation
fig,[ax1,ax2]=plt.subplots(1,2,figsize=(8,10))
times=np.linspace(0,1000-240,26)
alt=ftcp_c['altitude'][0:150]/1000

def animate(i):
    ax1.clear()
    ax2.clear()

    ax1.plot(rtm_hi[i,:],alt,c='r')
    ax1.plot(rtm_lo[i,:],alt,c='b')
    ax1.set_title("rtm")
    #ax1.set_xlim(300,330)

    ax2.plot(T_hi[i,:],alt,c='r')
    ax2.plot(T_lo[i,:],alt,c='b')
    ax2.set_title("T_in_K")
    #ax2.set_xlim(260,310)


    return fig
ani=FuncAnimation(fig,animate,frames=26,interval=100,repeat=True)
HTML(ani.to_jshtml())

# %%
from matplotlib import animation
fig,[ax1,ax2,ax3]=plt.subplots(1,3,figsize=(12,10))
times=np.linspace(0,950,100)
alt=fpzt['altitude'][:]/1000
mean=fpzt['thvm'][times.astype(int),:,0,0]
data_hi=fpzthi['thvm'][times.astype(int),:,0,0]#-mean
data_lo=fpztlo['thvm'][times.astype(int),:,0,0]#-mean
grad_hi=np.gradient(data_hi,axis=1)
grad_lo=np.gradient(data_lo,axis=1)
grad2_hi=np.gradient(grad_hi,axis=1)
grad2_lo=np.gradient(grad_lo,axis=1)
T_hi=fpzthi['T_in_K'][times.astype(int),:,0,0]
T_lo=fpztlo['T_in_K'][times.astype(int),:,0,0]
rtm_hi=fpzthi['rtm'][times.astype(int),:,0,0]
rtm_lo=fpztlo['rtm'][times.astype(int),:,0,0]
def animate(i):
    ax1.clear()
    ax2.clear()
    ax3.clear()
    
    ax1.plot(data_hi[i,:],alt,c='r')
    ax1.plot(data_lo[i,:],alt,c='b')
    ax1.set_title("VPT")
    ax1.set_xlim(300,330)
    ax1.set_ylim(0,7)
    
    ax2.plot(T_hi[i,:],alt,c='r')
    ax2.plot(T_lo[i,:],alt,c='b')
    ax2.set_title("T_in_K")
    ax2.set_xlim(260,310)
    ax2.set_ylim(0,7)
    
    ax3.plot(rtm_hi[i,:],alt,c='r')
    ax3.plot(rtm_lo[i,:],alt,c='b')
    ax3.set_title("rtm")
    #ax3.set_xlim(-.1,.4)
    ax3.set_ylim(0,7)
    
    return fig
    

ani=FuncAnimation(fig,animate,frames=100,interval=100,repeat=True)
#FFwriter = animation.FFMpegWriter(fps=10)
#ani.save('ani_1.mp4',writer=FFwriter)
from IPython.display import HTML
HTML(ani.to_jshtml())


# %%
def calc_maxc(vpt_h,cc1):
    maxc=np.zeros((vpt_h.shape[0],))
    for t in range(vpt_h.shape[0]):
        vptg =np.gradient(vpt_h[t,:],)
        vptg[0:5]=0
        vptg=np.cumsum(vptg)
        maxc[t]=np.argmin(np.abs(vptg-cc1))/150*6
    return maxc
def calc_maxg(vpt_h,cc1):
    maxc=np.zeros((vpt_h.shape[0],))
    for t in range(vpt_h.shape[0]):
        vptg =np.gradient(vpt_h[t,:],)
        vptg[0:5]=0
        maxc[t]=np.argmax(vptg)/150*6
    return maxc
def calc_blh(rig,lim):
    maxc=np.zeros((rig.shape[0],))
    for t in range(vpt_h.shape[0]):
        ri_=rig[t,:]
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
        maxc[t]=idx/150*6
    return maxc
    


# %%
plt.plot(fs1c['lh'][:,0,0,0])
plt.plot(fs2c['lh'][:,0,0,0])
plt.plot(fscp['lh'][:,0,0,0])
plt.plot(fs2c_c['lh'][:,0,0,0],'k-o')

# %%
fc1c=nc.Dataset(dir1c+'k_1/clusters.nc')
fc2c=nc.Dataset(dir2c+'k_2/clusters.nc')
fccp=nc.Dataset(dircp+'k_2/clusters.nc')

# %%
#plt.imshow(fc1c['cluster'][10,:])
plt.imshow(fc2c['cluster'][10,:])

# %%
fpsfc=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')

# %%
for v in fc1c.variables:
    print(v)

# %%

# %%

# %%

# %%

# %%
frc_c=read_forcings(tdir+'k_2/c_1/arm_forcings_f.in')
frc_h=read_forcings(tdir+'k_2/c_2/arm_forcings_f.in')

# %%
frc_o=read_forcings(ccdir+'k_2/c_1/arm_forcings_f.in')

# %%
frc_o.keys()


# %%

# %%

# %%

# %%
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


# %%

# %%
day = 20180707
dir2c = clubb_dir+'sgp_2c_p/sgp_'+str(day)+'/'
dir1c = clubb_dir+'sgp_1c_p/sgp_'+str(day)+'/'
dircp = clubb_dir+'sgp_cpl_p/sgp_'+str(day)+'/'
frc_cc= read_forcings(dircp+'k_2/c_1/arm_forcings_f.in')
frc_ch= read_forcings(dircp+'k_2/c_2/arm_forcings_f.in')
frc_0 = read_forcings(dir2c+'k_2/c_1/arm_forcings_f.in')

ftcp_c=nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
ftcp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
ft2c_c=nc.Dataset(dir2c+'k_2/c_1/output/arm_zt.nc','r')
ft2c_h=nc.Dataset(dir2c+'k_2/c_2/output/arm_zt.nc','r')
ft1c=nc.Dataset(dir1c+'k_1/c_1/output/arm_zt.nc','r')

fig=plt.figure(figsize=(10,5))

ax=plt.subplot(2,2,1)
tfrc=frc_ch['T_f[K\s]'][24:-24,0:150]+frc_cc['T_f[K\s]'][24:-24:,0:150]-2*frc_0['T_f[K\s]'][24:-24:,0:150]
im=plt.imshow(tfrc.T,origin='lower',cmap='coolwarm',vmin=-np.max(np.abs(tfrc)),vmax=np.max(np.abs(tfrc)),extent=(7,20,0,6))
plt.title('T_forcing -- '+str(day))
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im,cax=cax)

ax=plt.subplot(2,2,2)
tfrc=frc_ch['rtm_f[kg\kg\s]'][24:-24,0:150]+frc_cc['rtm_f[kg\kg\s]'][24:-24:,0:150]-2*frc_0['rtm_f[kg\kg\s]'][24:-24:,0:150]
im=plt.imshow(tfrc.T,origin='lower',cmap='coolwarm',vmin=-np.max(np.abs(tfrc)),vmax=np.max(np.abs(tfrc)),extent=(7,20,0,6))
plt.title('RTM_forcing -- '+str(day))
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im,cax=cax)

ax=plt.subplot(2,2,4)
data=ftcp_h['rtm_ma'][120:-120,0:150,0,0]-ft1c['rtm_ma'][120:-120,0:150,0,0]
im=plt.imshow(data.T,origin='lower',cmap='coolwarm',vmin=-np.max(np.abs(data)),vmax=np.max(np.abs(data)),extent=(7,20,0,6))
plt.title('RTM_forcing -- '+str(day))
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im,cax=cax)

ax=plt.subplot(2,2,3)
data=ftcp_h['thlm_ma'][120:-120,0:150,0,0]-ft1c['thlm_ma'][120:-120,0:150,0,0]
im=plt.imshow(data.T,origin='lower',cmap='coolwarm',vmin=-np.max(np.abs(data)),vmax=np.max(np.abs(data)),extent=(7,20,0,6))
plt.title('THLM_forcing -- '+str(day))
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im,cax=cax)

# %%

# %%

# %%
fig=plt.figure(figsize=(20,5))
im=plt.imshow(ftcp_h['thlm_forcing'][600:700,0:10,0,0].T,origin='lower')
ax=plt.gca()
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im,cax=cax)

# %%
fig=plt.figure(figsize=(20,5))
im=plt.imshow(ftcp_h['rtm_ma'][600:700,0:50,0,0].T,origin='lower')
ax=plt.gca()
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im,cax=cax)

# %%
fig=plt.figure(figsize=(20,5))
im=plt.imshow(ftcp_h['wm'][600:700,0:50,0,0].T,origin='lower')
ax=plt.gca()
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im,cax=cax)

# %%
np.mean(ftcp_h['thlm'][600:700,0:50,0,0])

# %%
plt.imshow(fmcp_h['thlp2_ma'][:,:,0,0].T,origin='lower')
plt.colorbar()

# %%
for k in ftcp_h.variables:
    print(k)

# %%
h2=[]
h3=[] #=h6
h4=[]
h5=[]
h2=[] #=h7
hr=15
dT=3
for day in days:
    dir2c = clubb_dir+'sgp_2c_p/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_p/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_p/sgp_'+str(day)+'/'
    
    ftcp_c=nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
    ftcp_h=nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
    alt=ftcp_c['altitude'][0:150]
    vptc = ftcp_c['thvm'][(hr-5)*60,0:150,0,0]
    vpth = ftcp_h['thvm'][(hr-5)*60,0:150,0,0]
    h4.append(np.argmax((vpth-(np.min(vpth)+dT))>0))
    if vpth[h4[-1]]>=vptc[h4[-1]]:
        h4[-1]=0
        h5.append(0)
        h3.append(0)
        h2.append(0)
    else:
        h5.append(np.argmax((vptc-vpth[h4[-1]])>0)-1)
        h3.append(np.argmax((vpth-vptc)<0)-1)
        h2.append(min(h5[-1]-h3[-1],h3[-1]-1))
    '''
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(vpth-vptc,alt,'b-')
    #plt.plot(vpth,alt,'r-')
    plt.title(day)
    plt.subplot(1,2,2)
    plt.plot([0,1],np.array([h2[-1],h2[-1]])*40)
    plt.plot([0,1],np.array([h3[-1],h3[-1]])*40)
    plt.plot([0,1],np.array([h4[-1],h5[-1]])*40)
    plt.ylim(0,6000)
    '''
    print('.',end='')
dta=(np.array(h4)-np.array(h5))
plt.hist(dta*40)
print()

# %%

# %%
import pickle
fp=open('test.p','rb')
data = pickle.load(fp)
data2 = pickle.load(fp)
fp.close()

# %%
data2

# %%
fp.close()

# %%
