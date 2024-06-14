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
import scipy as sci
import rasterio
import scipy.ndimage
sns.set_theme()
plt.rcParams.update({'figure.max_open_warning': 0})

# %%
clubb_dir='/home/tsw35/tyche/clubb/'

# %%
day='20170717'
dirr=clubb_dir+'mar27fit_cpl/sgp_'+day+'/'
#dirr=clubb_dir+'test_cpl/'
ft=nc.Dataset(dirr+'k_2/agg_outzt.nc','r')
fs=nc.Dataset(dirr+'k_2/agg_outsfc.nc','r')
fm=nc.Dataset(dirr+'k_2/agg_outzm.nc','r')
ftc=nc.Dataset(dirr+'k_2/c_1/output/arm_zt.nc','r')
fsc=nc.Dataset(dirr+'k_2/c_1/output/arm_sfc.nc','r')
fmc=nc.Dataset(dirr+'k_2/c_1/output/arm_zm.nc','r')
fth=nc.Dataset(dirr+'k_2/c_2/output/arm_zt.nc','r')
fsh=nc.Dataset(dirr+'k_2/c_2/output/arm_sfc.nc','r')
fmh=nc.Dataset(dirr+'k_2/c_2/output/arm_zm.nc','r')
fc=nc.Dataset(dirr+'k_2/clusters.nc','r')

dirr1=clubb_dir+'mar27fit_2c/sgp_'+day+'/'
#dirr1=clubb_dir+'tall_1c/sgp_20160625/'
ft1=nc.Dataset(dirr1+'k_2/agg_outzt.nc','r')
fs1=nc.Dataset(dirr1+'k_2/agg_outsfc.nc','r')
fm1=nc.Dataset(dirr1+'k_2/agg_outzm.nc','r')

#dirr1=clubb_dir+'mar27fit_1c/sgp_'+day+'/'
#dirr1=clubb_dir+'tall_1c/sgp_20160625/'
#ft1=nc.Dataset(dirr1+'k_1/c_1/output/arm_zt.nc','r')
#fs1=nc.Dataset(dirr1+'k_1/c_1/output/arm_sfc.nc','r')
#fm1=nc.Dataset(dirr1+'k_1/c_1/output/arm_zm.nc','r')

# %%
ftv=['thlm','um','rtm','thvm','wm','cloud_cover','rtm_forcing','thlm_forcing']
fsv=['sh','lh','lwp']
fmv=['wp2','up2','vp2','Richardson_num','bv_freq_sqd','wprtp','wpthlp','thlp2','rtp2']

# %%
fig,axs=plt.subplots(4,6,figsize=(24,16),dpi=200)
i=0
j=0
for v in ftv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=ft[v][:,0:,0,0]
    im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
    fig.colorbar(im,ax=ax)
    ax.set_title(v)
for v in fmv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=fm[v][:,0:,0,0]
    im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
    fig.colorbar(im,ax=ax)
    ax.set_title(v)
for v in fsv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=fs[v][:,0,0,0]
    ax.plot(data)
    ax.set_title(v)
ax=axs[i,j]
j=j+1

data=fc['wbar'][0,:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
plt.title('wbar_c')

ax=axs[i,j]
j=j+1
data=fc['wbar'][1,:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
plt.title('wbar_h')

ax=axs[i,j]
j=j+1
data=fc['u_r'][:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
plt.title('u_r')

ax=axs[i,j]
j=j+1
data=fc['z_circh'][:]
ax.plot(data)
plt.title('z_circh')

# %%

# %%

# %%

# %%

# %%

# %%
fig,axs=plt.subplots(4,6,figsize=(24,16),dpi=200)
i=0
j=0
for v in ftv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=ft[v][:,0:150,0,0]-ft1[v][:,0:150,0,0]
    abmax=np.percentile(np.abs(data),99)
    im=ax.imshow(data.T,vmin=-abmax,vmax=abmax,origin='lower',cmap='coolwarm',extent=(5,22,0,12))
    fig.colorbar(im,ax=ax)
    ax.set_title(v)
for v in fmv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=fm[v][:,0:150,0,0]-fm1[v][:,0:150,0,0]
    abmax=np.percentile(np.abs(data),99)
    im=ax.imshow(data.T,vmin=-abmax,vmax=abmax,origin='lower',cmap='coolwarm',extent=(5,22,0,12))
    fig.colorbar(im,ax=ax)
    ax.set_title(v)
for v in fsv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=fs[v][:,0,0,0]-fs1[v][:,0,0,0]
    ax.plot(data)
    ax.set_title(v)
ax=axs[i,j]
j=j+1
data=fc['wbar'][0,:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
plt.title('wbar_c')

ax=axs[i,j]
j=j+1
data=fc['wbar'][1,:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
ax.set_title('wbar_h')

ax=axs[i,j]
j=j+1
data=fc['u_r'][:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
ax.set_title('u_r')

ax=axs[i,j]
j=j+1
data=fc['z_circh'][:]
ax.plot(data)
ax.set_title('z_circh')

# %%
fig,axs=plt.subplots(4,6,figsize=(24,16),dpi=200)
i=0
j=0
for v in ftv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=fth[v][:,0:150,0,0]-ftc[v][:,0:150,0,0]
    abmax=np.percentile(np.abs(data),99)
    im=ax.imshow(data.T,vmin=-abmax,vmax=abmax,origin='lower',cmap='coolwarm',extent=(5,22,0,12))
    fig.colorbar(im,ax=ax)
    ax.set_title(v)
for v in fmv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=fmh[v][:,0:150,0,0]-fmc[v][:,0:150,0,0]
    abmax=np.percentile(np.abs(data),99)
    im=ax.imshow(data.T,vmin=-abmax,vmax=abmax,origin='lower',cmap='coolwarm',extent=(5,22,0,12))
    fig.colorbar(im,ax=ax)
    ax.set_title(v)
for v in fsv:
    ax=axs[i,j]
    if j<5:
        j=j+1
    else:
        i=i+1
        j=0
    data=fsh[v][:,0,0,0]-fsc[v][:,0,0,0]
    ax.plot(data)
    ax.set_title(v)

ax=axs[i,j]
j=j+1
data=fc['wbar'][0,:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
plt.title('wbar_c')

ax=axs[i,j]
j=j+1
data=fc['wbar'][1,:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
plt.title('wbar_h')

ax=axs[i,j]
j=j+1
data=fc['u_r'][:,:150]
im=ax.imshow(data.T,origin='lower',cmap='terrain',extent=(5,22,0,12))
fig.colorbar(im,ax=ax)
plt.title('u_r')

ax=axs[i,j]
j=j+1
data=fc['z_circh'][:]
ax.plot(data)
plt.title('z_circh')


# %%
import pickle
fp=open('tunefeb23.p','rb')
mean_val=pickle.load(fp)
peak_val=pickle.load(fp)
h_c_lwpm=pickle.load(fp)
mean_rtm=pickle.load(fp)
h_c_rtm=pickle.load(fp)
crs=pickle.load(fp)
incs=pickle.load(fp)
cc1s=pickle.load(fp)
lwp =pickle.load(fp)
rtp2 = pickle.load(fp)
maxzl = pickle.load(fp)
maxzh = pickle.load(fp)
urmu = pickle.load(fp)
urmx = pickle.load(fp)
urfc = pickle.load(fp)
ccc0 = pickle.load(fp)
thlp2 = pickle.load(fp)
wp2 = pickle.load(fp)
met = pickle.load(fp)
fp.close()

# %%
days=list(mean_val.keys())
days.sort()


# %%
def plot_test(param,test,name,c2=0,day=0):
    param=np.copy(param)
    #param=param[cc1s[day]==3]
    #test=test[cc1s[day]==3]
    gu=np.unique(param)
    gu.sort()
    test_line=[]
    for g in gu:
        test_line.append(np.nanmean(test[param==g]))
    try:
        c2 = c2+1
        plt.scatter(param,test,alpha=.25)
    except:
        c2=np.array(c2)
        c2=c2[cc1s[day]==3]
        color=plt.cm.coolwarm((np.array(c2)-.5)/(1.5-.5))
        plt.scatter(param,test,c=color,alpha=.25)
    #plt.plot(gu,test_line,'b-')
    #xaxh=np.nanmax(gu[gu>0])-np.nanmin(gu[gu>0])
    #plt.xlim(np.nanmin(gu[gu>0])-xaxh*.1,np.nanmax(gu[gu>0])+xaxh*.3)
    yaxh=np.nanmax(test[test>0])-np.nanmin(test[test>0])
    ymin=np.nanmin(test[test>0])-.1*yaxh
    ymax=np.nanmax(test[test>0])+.3*yaxh
    plt.ylim(ymin,ymax)
    plt.title(name)


# %%
import warnings

# %%
days=list(mean_val.keys())
days.sort()
import warnings
warnings.filterwarnings("ignore")
for day in days:
    
    plt.figure(figsize=(24,12))
    
    # CRS
    plt.subplot(3,8,1)
    plot_test(crs[day],wp2[day],'WP2')
    plt.ylabel('CR')
    plt.subplot(3,8,2)
    plot_test(crs[day],mean_val[day],'Mean LWP')
    plt.subplot(3,8,3)
    plot_test(crs[day],h_c_lwpm[day],'H/C Ratio LWP')
    plt.subplot(3,8,4)
    plot_test(crs[day],rtp2[day],'RTP2')
    plt.subplot(3,8,5)
    plot_test(crs[day],thlp2[day],'THLP2')
    plt.subplot(3,8,6)
    plot_test(crs[day],peak_val[day],'PEAK LWP')
    plt.subplot(3,8,7)
    #plot_test(crs[day],urmx[day],'MAX UR')
    
    plt.subplot(3,8,8)
    colors=['red','white','gainsboro','lightgrey','lightgrey','grey','darkgrey','dimgrey','black','black','black']
    gs=np.unique(crs[day])
    gs.sort() 
    leglist=[]
    for idx in range(len(gs)):
        g=gs[idx]
        c=colors[idx]
        m_LWP=np.mean(lwp[day][crs[day]==g,:],axis=0)
        plt.plot(m_LWP,color=c)
        leglist.append(str(g))
    plt.xlim(120,1000)
    plt.legend(leglist,loc='upper left',labelspacing=.25,borderaxespad=0.1)
    plt.title('Mean LWP in Time')
  
    
    plt.suptitle(day)
    plt.subplots_adjust(wspace=.2,hspace=.35)
    
    # C1S
    plt.subplot(3,8,9)
    plot_test(cc1s[day],wp2[day],'WP2')
    plt.ylabel('CR')
    plt.subplot(3,8,10)
    plot_test(cc1s[day],mean_val[day],'Mean LWP')
    plt.subplot(3,8,11)
    plot_test(cc1s[day],h_c_lwpm[day],'H/C Ratio LWP')
    plt.subplot(3,8,12)
    plot_test(cc1s[day],rtp2[day],'RTP2')
    plt.subplot(3,8,13)
    plot_test(cc1s[day],thlp2[day],'THLP2')
    plt.subplot(3,8,14)
    plot_test(cc1s[day],peak_val[day],'PEAK LWP')
    plt.subplot(3,8,15)
    #plot_test(cc1s[day],urmx[day],'MAX UR')
    
    plt.subplot(3,8,16)
    colors=['red','white','gainsboro','lightgrey','lightgrey','grey','darkgrey','dimgrey','black','black','black']
    gs=np.unique(cc1s[day])
    gs.sort() 
    leglist=[]
    for idx in range(len(gs)):
        g=gs[idx]
        c=colors[idx]
        m_LWP=np.mean(lwp[day][cc1s[day]==g,:],axis=0)
        plt.plot(m_LWP,color=c)
        leglist.append(str(g))
    plt.xlim(120,1000)
    plt.legend(leglist,loc='upper left',labelspacing=.25,borderaxespad=0.1)
    plt.title('Mean LWP in Time')
  
    
    plt.suptitle(day)
    plt.subplots_adjust(wspace=.2,hspace=.35)
    
    # MET
    met2=met[day]
    met3=np.zeros(met2.shape)
    met3[met2=='denv_a']=0
    met3[met2=='denv_s']=1
    met3[met2=='denv_ds']=2
    plt.subplot(3,8,17)
    plot_test(met3,wp2[day],'WP2')
    plt.ylabel('CR')
    plt.subplot(3,8,18)
    plot_test(met3,mean_val[day],'Mean LWP')
    plt.subplot(3,8,19)
    plot_test(met3,h_c_lwpm[day],'H/C Ratio LWP')
    plt.subplot(3,8,20)
    plot_test(met3,rtp2[day],'RTP2')
    plt.subplot(3,8,21)
    plot_test(met3,thlp2[day],'THLP2')
    plt.subplot(3,8,22)
    plot_test(met3,peak_val[day],'PEAK LWP')
    plt.subplot(3,8,23)
    #plot_test(met3,urmx[day],'MAX UR')
    
    plt.subplot(3,8,24)
    colors=['red','white','gainsboro','lightgrey','lightgrey','grey','darkgrey','dimgrey','black','black','black']
    gs=np.unique(met3)
    gs.sort() 
    leglist=[]
    for idx in range(len(gs)):
        g=gs[idx]
        c=colors[idx]
        m_LWP=np.mean(lwp[day][met3==g,:],axis=0)
        plt.plot(m_LWP,color=c)
        leglist.append(str(g))
    plt.xlim(120,1000)
    plt.legend(leglist,loc='upper left',labelspacing=.25,borderaxespad=0.1)
    plt.title('Mean LWP in Time')
  
    
    plt.suptitle(day)
    plt.subplots_adjust(wspace=.2,hspace=.35)

# %%

# %%
testdir=clubb_dir+'mar15fit_cpl/'
dirles='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'

for file in filelist:
    day=file[4:]
    ftc=nc.Dataset(clubb_dir+'mar15fit_cpl/'+file+'/k_2/c_1/output/arm_zt.nc','r')
    fth=nc.Dataset(clubb_dir+'mar15fit_cpl/'+file+'/k_2/c_2/output/arm_zt.nc','r')
    plt.figure(figsize=(3,3))
    data=fth['thvm'][:,0:50,0,0]-ftc['thvm'][:,0:50,0,0]
    plt.imshow(np.abs(data.T),origin='lower',cmap='coolwarm',extent=(5,22,0,20),vmin=-1,vmax=1)
    plt.colorbar()
    plt.title(day)
    break

# %%
day='20160625'
fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
data=fplm['u'][:,166:]
plt.imshow(np.abs(data.T),origin='lower',cmap='terrain',extent=(7,22,0,6),vmin=-4,vmax=10)
plt.colorbar()

# %%
plt.plot(fplm['z'][1,166:])

# %%
testdir=clubb_dir+'mar20fit_cpl/'
dirles='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'

for file in filelist:
    day=file[4:]
    ft1=nc.Dataset(clubb_dir+'mar20fit_1c/'+file+'/k_1/c_1/output/arm_zt.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    
    #data=np.sqrt(ft1['um'][120:,0:125,0,0]**2+ft1['vm'][120:,0:125,0,0]**2)
    #data3=np.sqrt(fplm['u'][:,0:166]**2+fplm['v'][:,0:166]**2)
    
    data=ft1['um'][:,0:125,0,0]
    data2=ft1['vm'][:,0:125,0,0]
    data3=fplm['u'][:,0:166]
    data4=fplm['v'][:,0:166]
    
    vmax=max(np.max(np.abs(data)),np.max(np.abs(data3)))
    vmax=max(vmax,np.max(np.abs(data2)),np.max(np.abs(data4)))
    
    plt.figure(figsize=(6,6))
    plt.subplot(2,2,1)
    im=plt.imshow(np.abs(data.T),origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    plt.colorbar()
    plt.title('CLUBB U: '+day)
    
    plt.subplot(2,2,2)
    
    plt.imshow(np.abs(data3.T),origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    plt.colorbar()
    plt.title('LES U: '+day)
    
    plt.subplot(2,2,3)
    im=plt.imshow(np.abs(data2.T),origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    plt.colorbar()
    plt.title('CLUBB V: '+day)
    
    plt.subplot(2,2,4)
    
    plt.imshow(np.abs(data4.T),origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    plt.colorbar()
    plt.title('LES V: '+day)
    

# %%
import re

# %%
urr=fc['u_rr'][:]

# %%
plt.imshow(urr[1,:,:].T,origin='lower',cmap='terrain',extent=(7,23,0,10),vmin=0,vmax=1)
plt.colorbar()


# %%
def true_u(ur,ur0,uu,vv,wfrac):
    vv2=np.zeros((192,125))
    uu2=np.zeros((192,125))
    for t in range(0,955,5):
        vv2[int(t/5),:]=vv[t,:]
        uu2[int(t/5),:]=uu[t,:]
    vmsk=vv2<ur0
    umsk=uu2<ur0
    urout=ur.copy()
    urout[umsk &~vmsk]=urout[umsk &~vmsk]/wfrac[0]
    urout[vmsk &~umsk]=urout[vmsk &~umsk]/wfrac[1]
    return urout


# %%
ur.shape


# %%
def true_u2(ur,urr,wfrac):
    urout=ur.copy()
    for t in range(ur.shape[0]):
        urout[(urr[0,t,:]>=.99)]=urout[(urr[0,t,:]>=.99)]/wfrac[1]


# %%
testdir=clubb_dir+'apr10_localfit_cpl/'
dirles='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'

for file in filelist:
    vmax=4
    day=file[4:]
    fs1=nc.Dataset(clubb_dir+'apr10_localfit_1c/'+file+'/k_1/c_1/output/arm_sfc.nc','r')
    fs2=nc.Dataset(clubb_dir+'apr10_localfit_2c/'+file+'/k_2/agg_outsfc.nc','r')
    fsc=nc.Dataset(testdir+file+'/k_2/agg_outsfc.nc','r')
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    print(np.max(fc['wbar'][:]))

# %%
plt.figure(figsize=(15,8))
warnings.filterwarnings("ignore")
leglist=[]
for file in os.listdir(clubb_dir):
    if 'apr13' not in file:
        continue
    if 'cpl' not in file:
        continue
    if '10' in file:
        continue
    if '20' in file:
        continue
    if '100' in file:
        subdir=file[0:15]
    else:
        subdir=file[0:14]
    #'apr13_dz60dt01'
    leglist.append(subdir[6:])
    testdir=clubb_dir+subdir+'_cpl/'
    dirles='/home/tsw35/soteria/clubb/data/les_param/'
    filelist=os.listdir(testdir)
    filelist.sort()
    les_1c ='/home/tsw35/tyche/data/LES_1C/'

    km5=83
    
    i=1
    for file in filelist:
        vmax=4
        day=file[4:]
        fs1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_sfc.nc','r')
        fs2=nc.Dataset(clubb_dir+subdir+'_2c/'+file+'/k_2/agg_outsfc.nc','r')
        fsc=nc.Dataset(testdir+file+'/k_2/agg_outsfc.nc','r')
        fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
        plt.subplot(2,5,i)
        plt.plot(fsc['lwp'][:,0,0,0],alpha=.5)
        #maxzh=np.max(fc['u_r'][:],axis=1)
        #maxzh[maxzh<=0]=float('nan')
        #plt.plot(maxzh,alpha=.5)
        plt.title(file)
        
        i=i+1
        
subdir='apr10_localfit'
leglist.append(subdir[6:])
testdir=clubb_dir+subdir+'_cpl/'
dirles='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'

km5=83

i=1
for file in filelist:
    vmax=4
    day=file[4:]
    fs1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_sfc.nc','r')
    fs2=nc.Dataset(clubb_dir+subdir+'_2c/'+file+'/k_2/agg_outsfc.nc','r')
    fsc=nc.Dataset(testdir+file+'/k_2/agg_outsfc.nc','r')
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    plt.subplot(2,5,i)
    plt.plot(fsc['lwp'][:,0,0,0],alpha=.5)
    #maxzh=np.max(fc['u_r'][:],axis=1)
    #maxzh[maxzh<=0]=float('nan')
    #plt.plot(maxzh,alpha=.5)
    plt.title(file)

    i=i+1
        
        
plt.subplot(2,5,10)
plt.legend(leglist)

# %%
alt=fplm['z'][1,:]
plt.plot(alt)
alt=fplm['z'][45,:]
plt.plot(alt)
alt=fplm['z'][90,:]
plt.plot(alt)
plt.show()

# %%
alt.shape

# %%
subdir='apr25_dz60dt06'
testdir=clubb_dir+subdir+'_cpl/'
dirles='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'

km5=83
km5l=166

for file in filelist:
    vmax=4
    day=file[4:]
    fs1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_sfc.nc','r')
    fs2=nc.Dataset(clubb_dir+subdir+'_2c/'+file+'/k_2/agg_outsfc.nc','r')
    fsc=nc.Dataset(testdir+file+'/k_2/agg_outsfc.nc','r')
    fscc=nc.Dataset(testdir+file+'/k_2/c_1/output/arm_sfc.nc','r')
    fsch=nc.Dataset(testdir+file+'/k_2/c_2/output/arm_sfc.nc','r')
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    ft=nc.Dataset(testdir+file+'/k_2/agg_outzt.nc','r')
    fm=nc.Dataset(testdir+file+'/k_2/agg_outzm.nc','r')
    fm1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_zm.nc','r')
    ft1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_zt.nc','r')
    
    fp=nc.Dataset(dirles+day+'.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    
    runpm=open(testdir+file+'/tw_run_param.txt','r').readlines()
    
    for ln in runpm:
        if 'c_ur' in ln:
            c_ur=re.findall('\d+\.\d+',ln)[0]
            
    ur0 =fc['u_r0'][:,:km5]
    uu  =ft['um'][:,:km5,0,0]
    vv  =ft['vm'][:,:km5,0,0]
    wfrac=fc['W_frac'][:]
    ur=fc['u_r'][:,:km5]
    ur[ur==-1]=float('nan')
    ur[ur==0]=float('nan')
    #ur0[ur0==-1]=float('nan')
    #ur0[ur0==0]=float('nan')
    maxzh=fc['z_circh'][:]
    
    data=ur#ur0#true_u(ur,ur0,uu,vv,wfrac)
    
    plt.figure(figsize=(12,14))
    plt.subplot(5,2,1)
    im=plt.imshow(np.abs(data.T),origin='lower',cmap='terrain',extent=(7,23,0,10),vmin=0,vmax=vmax)
    plt.colorbar()
    plt.plot(np.linspace(7,23,192),maxzh/1000*2,'k--')
    plt.title('CLUBB CIRC.: '+day+' : '+str(c_ur))
    
    
    
    
    plt.subplot(5,2,2)
    u90=fp['u90'][:]
    u10=fp['u10'][:]
    v90=fp['v90'][:]
    v10=fp['v10'][:]

    data1=u90[:,0:km5l].T-u10[:,0:km5l].T
    nch=(u90[:,0:km5l].T/np.abs(u90[:,0:km5l].T))==(u10[:,0:km5l].T/np.abs(u10[:,0:km5l].T))
    data4=(data1/2)**2
    data1[nch]=0
    data1=data1/2

    data2=v90[:,0:km5l].T-v10[:,0:km5l].T
    nch=(v90[:,0:km5l].T/np.abs(v90[:,0:km5l].T))==(v10[:,0:km5l].T/np.abs(v10[:,0:km5l].T))
    data4=data4+(data2/2)**2
    data2[nch]=0
    data2=data2/2
    
    '''
    data1=u90[:,0:km5l].T
    data1[np.abs(data1)<np.abs(u10[:,0:km5l].T)]=u10[:,0:km5l].T[np.abs(data1)<np.abs(u10[:,0:km5l].T)]
    
    data2=v90[:,0:km5l].T
    data2[np.abs(data2)<np.abs(v10[:,0:km5l].T)]=v10[:,0:km5l].T[np.abs(data2)<np.abs(v10[:,0:km5l].T)]
    '''
    
    data3=np.sqrt(data1**2+data2**2)
    data4=np.sqrt(data4)
    data3[data3==0]=float('nan')
    data4[data4<=.1]=float('nan')
    
    plt.imshow(data4,alpha=.2,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    plt.imshow(data3,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    
    plt.colorbar()
    plt.title('LES CIRC.: '+day)

    
    
    plt.subplot(5,2,3)
    lwpt=fs1['lwp'][:,0,0,0]
    lwpm=fs2['lwp'][:,0,0,0]
    lwpc=fsc['lwp'][:,0,0,0]
    plt.plot(np.linspace(7,23,len(lwpm)),lwpm,'b--')
    plt.plot(np.linspace(7,23,len(lwpm)),lwpt,'r--')
    plt.plot(np.linspace(7,23,len(lwpm)),lwpc,'r-')
    plt.title('CLUBB LWP')
    
    
    
    plt.subplot(5,2,4)
    lwpt=fplt['lwp'][:]
    lwpm=fplm['lwp'][:]
    plt.plot(np.linspace(7,22,len(lwpm)),lwpm,'b--')
    plt.plot(np.linspace(7,22,len(lwpm)),lwpt,'r-')
    plt.title('LES LWP')
    
    
    data1=fm['thlp2'][:,:km5,0,0].T#-fm1['thlp2'][:,:km5,0,0].T
    data1=sci.ndimage.filters.gaussian_filter(data1,sigma=5)
    data2=fplt['thl2'][:,:km5l].T#-fplm['thl2'][:,:km5l].T
    vmax=max(np.max(data1),np.max(data2[15:,:]))
    vmin=0
    
    plt.subplot(5,2,5)
    plt.imshow(data1,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=vmin,vmax=vmax)
    plt.colorbar()
    
    plt.subplot(5,2,6)
    plt.imshow(data2,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=vmin,vmax=vmax)
    plt.colorbar()
    
    
    data1m=fscc['lwp'][:,0,0,0]
    data1t=fsch['lwp'][:,0,0,0]
    data2t=np.sum(.5*(fplt['u2'][:,:km5l]+fplt['u2'][:,:km5l]+fplt['u2'][:,:km5l]),axis=1)
    data2m=np.sum(.5*(fplm['u2'][:,:km5l]+fplm['u2'][:,:km5l]+fplm['u2'][:,:km5l]),axis=1)
    vmax=max(np.max(data1),np.max(data2[15:,:]))
    vmin=0 #-vmax
    
    plt.subplot(5,2,7)
    plt.plot(data1m)
    plt.plot(data1t)
    
    plt.subplot(5,2,8)
    plt.plot(data2m)
    plt.plot(data2t)
    
    
    data1=ft['rcm'][:,:km5,0,0].T
    data2=fplt['lwc'][:,:km5l].T
    vmax=max(np.percentile(data1,99),np.percentile(data2[:,:],99))
    vmin=0
    data1[data1==0]=float('nan')
    data2[data2==0]=float('nan')
    
    plt.subplot(5,2,9)
    plt.imshow(data1,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=vmin,vmax=vmax)
    plt.colorbar()
    
    plt.subplot(5,2,10)
    plt.imshow(data2,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=vmin,vmax=vmax)
    plt.colorbar()

# %%
print(filelist)

# %%
subdir='apr13_dz60dt06'
testdir=clubb_dir+subdir+'_cpl/'
dirles='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'

km5=83
km5l=166

for file in filelist:
    vmax=4
    day=file[4:]
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    ft=nc.Dataset(testdir+file+'/k_2/agg_outzt.nc','r')
    
    fp=nc.Dataset(dirles+day+'.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    
    runpm=open(testdir+file+'/tw_run_param.txt','r').readlines()
    
    for ln in runpm:
        if 'c_ur' in ln:
            c_ur=re.findall('\d+\.\d+',ln)[0]
            
    ur0 =fc['u_r0'][:,:km5]
    uu  =ft['um'][:,:km5,0,0]
    vv  =ft['vm'][:,:km5,0,0]
    wfrac=fc['W_frac'][:]
    ur=fc['u_r'][:,:km5]
    ur[ur==-1]=float('nan')
    ur[ur==0]=float('nan')
    #ur0[ur0==-1]=float('nan')
    #ur0[ur0==0]=float('nan')
    maxzh=fc['z_circh'][:]
    
    data=ur#ur0#true_u(ur,ur0,uu,vv,wfrac)
    
    plt.figure(figsize=(12,14))
    plt.subplot(2,2,1)
    im=plt.imshow(np.abs(data.T),origin='lower',cmap='terrain',extent=(7,23,0,10),vmin=0,vmax=vmax)
    plt.colorbar()
    plt.plot(np.linspace(7,23,192),maxzh/1000*2,'k--')
    plt.title('CLUBB CIRC.: '+day+' : '+str(c_ur))
    
    plt.subplot(2,2,2)
    u90=fp['u90'][:]
    u10=fp['u10'][:]
    v90=fp['v90'][:]
    v10=fp['v10'][:]

    data1=u90[:,0:km5l].T-u10[:,0:km5l].T
    nch=(u90[:,0:km5l].T/np.abs(u90[:,0:km5l].T))==(u10[:,0:km5l].T/np.abs(u10[:,0:km5l].T))
    data4=(data1/2)**2
    data1[nch]=0
    data1=data1/2

    data2=v90[:,0:km5l].T-v10[:,0:km5l].T
    nch=(v90[:,0:km5l].T/np.abs(v90[:,0:km5l].T))==(v10[:,0:km5l].T/np.abs(v10[:,0:km5l].T))
    data4=data4+(data2/2)**2
    data2[nch]=0
    data2=data2/2
    
    '''
    data1=u90[:,0:km5l].T
    data1[np.abs(data1)<np.abs(u10[:,0:km5l].T)]=u10[:,0:km5l].T[np.abs(data1)<np.abs(u10[:,0:km5l].T)]
    
    data2=v90[:,0:km5l].T
    data2[np.abs(data2)<np.abs(v10[:,0:km5l].T)]=v10[:,0:km5l].T[np.abs(data2)<np.abs(v10[:,0:km5l].T)]
    '''
    
    data3=np.sqrt(data1**2+data2**2)
    data4=np.sqrt(data4)
    data3[data3==0]=float('nan')
    data4[data4<=.1]=float('nan')
    
    plt.imshow(data4,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    #plt.imshow(data3,origin='lower',cmap='terrain',extent=(7,22,0,10),vmin=0,vmax=vmax)
    
    plt.colorbar()
    plt.title('LES CIRC.: '+day)


# %%
tskdir='/home/tsw35/soteria/clubb/data/surfaces_5k/lw'
tsk={}
sfc_dir   = '/home/tsw35/soteria/clubb/data/surfaces_5k/'
for file in filelist:
    file=file[4:]
    stdt=datetime.datetime(int(file[0:4]),int(file[4:6]),int(file[6:8]),12,0)
    lwg,lwv=read_sfc_data('lw',16,stdt)
    tskg=(lwg/(5.67*10**(-8)))**(1/4)
    tsk[file]=np.array([np.std(tskg,axis=(1,2)),np.std(tskg,axis=(1,2))])
    print('.',end='',flush=True)

# %%
dayst={}
for day in filelist:
    day=day[4:]
    dayst[day]=day[0:4]+'-'+day[4:6]+'-'+day[6:8]

# %%
subdir='apr13_dz60dt06'
testdir=clubb_dir+subdir+'_cpl/'
dir2c='/home/tsw35/soteria/clubb/data/les_param2/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'

km5=83
km5l=166
frc=1.25


cc1=.85
#filelist=['20160625.nc','20170716.nc','20170717.nc','20190707.nc','20180709.nc','20180707.nc','20170802.nc','20150801.nc','20160610.nc','20160716.nc']
for file in filelist:
    
    plt.figure(figsize=(5,4),dpi=300)
    
    vmax=4
    day=file[4:]
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    ft=nc.Dataset(testdir+file+'/k_2/agg_outzt.nc','r')
    
    fp=nc.Dataset(dirles+day+'.nc','r')

    runpm=open(testdir+file+'/tw_run_param.txt','r').readlines()
    
    for ln in runpm:
        if 'c_ur' in ln:
            c_ur=re.findall('\d+\.\d+',ln)[0]
            
    ur0 =fc['u_r0'][:,:km5]
    uu  =ft['um'][:,:km5,0,0]
    vv  =ft['vm'][:,:km5,0,0]
    wfrac=fc['W_frac'][:]
    ur=fc['u_r'][:,:km5]
    #ur[ur==-1]=0
    ur[ur==-1]=float('nan')
    ur[ur==0]=float('nan')
    ur=ur[0:-12,:]
    #ur0[ur0==-1]=float('nan')
    #ur0[ur0==0]=float('nan')
    maxzh=fc['z_circh'][:-12]
    data= np.abs(ur.T)
    if day=='20160625':
        data=data/fc['W_frac'][0]
    
    plt.subplot(2,1,1)
    
    plt.imshow(data,cmap='terrain',origin='lower',vmax=vmax,extent=(7,22,0,5*frc))
    #plt.colorbar(label=r'CPL $u_r$')
    #plt.title(dayst[day])
    plt.title('Model')
    
    plt.plot(np.linspace(7,22,len(maxzh)),maxzh/1000*frc,'k--')
    plt.ylim(0,frc*5)
    plt.xlim(8,20)
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([8,12,16,20],labels=[])
    plt.ylabel('Elevation ($km$)')
    #plt.xlabel('Local Time')
    
    
    
    
    
    plt.subplot(2,1,2)
    u90=fp['u90'][:]
    u10=fp['u10'][:]
    v90=fp['v90'][:]
    v10=fp['v10'][:]
    
    vpt=fp['vpt'][:]
    dtsk=2*tsk[day].T
    
    
    maxzh,maxzl,minzh = circ_h(vpt,dtsk,cc1)
    
    #plt.subplot(1,2,1)
    #plt.title(dayst[day])
    plt.title('LES')
    #data=u90[:,0:166].T-u10[:,0:166].T
    data=np.sqrt(((u90[:,0:166].T-u10[:,0:166].T)/2)**2+((v90[:,0:166].T-v10[:,0:166].T)/2)**2)
    
    nch1=(u90[:,0:km5l].T/np.abs(u90[:,0:km5l].T))==(u10[:,0:km5l].T/np.abs(u10[:,0:km5l].T))
    nch2=(v90[:,0:km5l].T/np.abs(v90[:,0:km5l].T))==(v10[:,0:km5l].T/np.abs(v10[:,0:km5l].T))
    
    
    data[nch1&nch2]=float('nan')
    
    if day=='20160625':
        data[int(3.5/5*km5l):,:]=float('nan')
        
    if day=='20170716':
        data[int(3/5*km5l):,:]=float('nan')
    
    #data=data/np.max(data,axis=0)
    
    plt.imshow(data,cmap='terrain',origin='lower',vmax=vmax,extent=(7,22,0,5*frc))
    #plt.colorbar(label=r'$u_r$ ($ms^{-1}$)')
    
    plt.plot(np.linspace(7,22,16),maxzh*30*frc/1000,'k--')
    #plt.plot(np.linspace(7,22,16),maxzl*30*frc/1000,'k-')
    #plt.plot(np.linspace(7,22,16),minzh*30*frc/1000,'k:')
    plt.ylim(0,frc*5)
    plt.xlim(8,20)
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    plt.ylabel('Elevation ($km$)')
    plt.xlabel('Local Time')
    plt.legend([r'$z_{max_1}$'],loc='upper left',framealpha=.9)

# %%
subdir='apr13_dz60dt06'
testdir=clubb_dir+subdir+'_cpl/'
dirles='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(testdir)
filelist.sort()
les_1c ='/home/tsw35/tyche/data/LES_1C/'
filelist=['sgp_20170717']

km5=83
km5l=166
frc=1.25
for file in filelist:
    vmax=4
    day=file[4:]
    fs1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_sfc.nc','r')
    fs2=nc.Dataset(clubb_dir+subdir+'_2c/'+file+'/k_2/agg_outsfc.nc','r')
    fsc=nc.Dataset(testdir+file+'/k_2/agg_outsfc.nc','r')
    fscc=nc.Dataset(testdir+file+'/k_2/c_1/output/arm_sfc.nc','r')
    fsch=nc.Dataset(testdir+file+'/k_2/c_2/output/arm_sfc.nc','r')
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    ft=nc.Dataset(testdir+file+'/k_2/agg_outzt.nc','r')
    fm=nc.Dataset(testdir+file+'/k_2/agg_outzm.nc','r')
    fm1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_zm.nc','r')
    ft1=nc.Dataset(clubb_dir+subdir+'_1c/'+file+'/k_1/c_1/output/arm_zt.nc','r')
    
    fp=nc.Dataset(dirles+day+'.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    
    runpm=open(testdir+file+'/tw_run_param.txt','r').readlines()
    
    for ln in runpm:
        if 'c_ur' in ln:
            c_ur=re.findall('\d+\.\d+',ln)[0]
            
    ur0 =fc['u_r0'][:,:km5]
    uu  =ft['um'][:,:km5,0,0]
    vv  =ft['vm'][:,:km5,0,0]
    wfrac=fc['W_frac'][:]
    ur=fc['u_r'][:,:km5]
    ur[ur==-1]=float('nan')
    ur[ur==0]=float('nan')
    #ur0[ur0==-1]=float('nan')
    #ur0[ur0==0]=float('nan')
    maxzh=fc['z_circh'][:]
    
    data=ur#ur0#true_u(ur,ur0,uu,vv,wfrac)
    
    plt.figure(figsize=(7,7),dpi=300)
    #plt.figure(figsize=(10,10),dpi=300)

    
    
    
    plt.subplot(3,2,1)
    lwpt=fs1['lwp'][:,0,0,0]
    lwpm=fs2['lwp'][:,0,0,0]
    lwpc=fsc['lwp'][:,0,0,0]
    #plt.plot(np.linspace(7,23,len(lwpm)),lwpm,'b--')
    plt.plot(np.linspace(7,23,len(lwpm)),lwpt,'b--')
    plt.plot(np.linspace(7,23,len(lwpm)),lwpc,'r-')
    plt.title('Model',fontsize=16)
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    plt.xlim(7,22)
    plt.ylim(0,0.042)
    plt.yticks([0,.01,.02,.03,.04])
    plt.legend(['Single Column','Coupled'])
    plt.ylabel('LWP ($kgm^{-2}$)')
    
    
    
    plt.subplot(3,2,2)
    lwpt=fplt['lwp'][:]
    lwpm=fplm['lwp'][:]
    plt.plot(np.linspace(7,22,len(lwpm)),lwpm,'b--')
    plt.plot(np.linspace(7,22,len(lwpm)),lwpt,'r-')
    plt.title('LES',fontsize=16)
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    plt.xlim(7,22)
    plt.ylim(0,0.042)
    plt.yticks([0,.01,.02,.03,.04],[])
    plt.legend(['Hmg','Het'])

    
    
    var1='thlp2'#'rtp2'
    var2='thl2'#'qv2'
    data1=fm[var1][:,:km5,0,0].T-fm1[var1][:,:km5,0,0].T
    data1=sci.ndimage.filters.gaussian_filter(data1,sigma=5)
    data2=fplt[var2][:,:km5l].T-fplm[var2][:,:km5l].T
    vmax=1.25
    vmin=-vmax
    
    cmap='coolwarm'
    plt.subplot(3,2,5)
    plt.imshow(data1,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
    #plt.colorbar()
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    plt.ylabel('Elevation ($km$)')
    
    plt.subplot(3,2,6)
    plt.imshow(data2,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
    plt.colorbar(label=r"$\Delta\theta_v^{\prime2}$ ($K^2$)")
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[])
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    
    
    
    cmap='Blues'
    data1=ft['rcm'][:,:km5,0,0].T*1000
    data2=fplt['lwc'][:,:km5l].T*1000
    vmax=max(np.percentile(data1,99),np.percentile(data2[:,:],99))
    vmin=0
    data1[data1<=.0003]=float('nan')
    data2[data2<=.0003]=float('nan')
    
    
    plt.subplot(3,2,3)
    plt.imshow(data1,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
    #plt.colorbar()
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[0,1,2,3,4,5])
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    plt.ylabel('Elevation ($km$)')
    
    plt.subplot(3,2,4)
    plt.imshow(data2,origin='lower',cmap=cmap,extent=(7,22,0,5*frc),vmin=vmin,vmax=vmax)
    plt.colorbar(label=r'Liquid Water ($g\ kg^{-1}$)')
    plt.yticks([0,frc,frc*2,frc*3,frc*4,frc*5],labels=[])
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    
    plt.subplots_adjust(hspace=.3)

# %%
fplt['lwc']

# %%
for file in filelist:
    
    plt.figure(figsize=(7,3),dpi=300)
    
    vmax=4
    day=file[4:]
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    plt.subplot(1,2,1)
    data=fc['u_rr'][0,:,0:50].T
    data[data<0]=float('nan')
    data[data>=.99]=float('nan')
    plt.imshow(data,origin='lower',cmap='terrain')
    plt.colorbar()
    plt.subplot(1,2,2)
    data=fc['u_rr'][1,:,0:50].T
    data[data<0]=float('nan')
    data[data>=.99]=float('nan')
    plt.imshow(data,origin='lower',cmap='terrain')
    plt.title(day)
    plt.colorbar()

# %%
fc['W_frac'][:]


# %%
def circ_h(vpt_,dtsk_,cc1=cc1):
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


# %%
fc.variables

# %%
fp0=nc.Dataset(testdir+'/sgp_20160625/k_2/c_1/output/arm_zt.nc','r')
fp1=nc.Dataset(testdir+'/sgp_20160625/k_2/c_2/output/arm_zt.nc','r')

plt.imshow(fp1['thvm'][:,0:125,0,0].T-fp0['thvm'][:,0:125,0,0].T,extent=(5,22,0,5),origin='lower',cmap='terrain',vmin=-2,vmax=2)
plt.colorbar()


# %%
for file in filelist:
    day=file[4:]
    fc=nc.Dataset(testdir+file+'/k_2/clusters.nc','r')
    data=fc['u_r0'][:,:125]-fc['u_r']
    data[data==-1]=float('nan')
    data[data==0]=float('nan')
    plt.imshow(np.abs(data.T),origin='lower',cmap='terrain',extent=(5,22,0,10),vmin=0,vmax=4)
    
    break

# %%
plt.imshow(np.abs(fp1['vm'][:,0:125,0,0].T),origin='lower',extent=(5,22,0,5),cmap='terrain')
plt.colorbar()

# %%
fp=nc.Dataset(dirles+'20160625.nc','r')
data=(fp['u90'][:]+fp['u10'][:])/2
plt.imshow(np.abs(data[:,0:166].T),origin='lower',extent=(7,22,0,5),cmap='terrain')
plt.colorbar()

# %%
fpp['thvm'][:].shape

# %%
for v in fc.variables:
    print(v)


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
