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
clubb_dir= '/home/tsw35/tyche/clubb/'

# %%
fp=nc.Dataset(clubb_dir+'test_cpl/k_2/clusters.nc','r')

# %%
fp.variables.keys()

# %%
data=fp['wbar'][0,:,:150]
print(data.shape)
data[data==0]=float('nan')
plt.imshow(data.T,vmin=-.4,vmax=.4,origin='lower',cmap='terrain',extent=(5,22,0,6))
plt.colorbar()

# %%
wbar=fp['wbar'][1,:,:150]
data=fp['u_r'][:,:150]
data[data==0]=float('nan')
#data[np.abs(data)<=0.1]=float('nan')
plt.imshow(data[:].T,vmin=-2,vmax=2,origin='lower',cmap='terrain',extent=(5,22,0,6))
plt.colorbar()

# %%
fs=nc.Dataset(clubb_dir+'test_cpl/k_2/agg_outsfc.nc','r')
ft=nc.Dataset(clubb_dir+'test_cpl/k_2/agg_outzt.nc','r')
ftc=nc.Dataset(clubb_dir+'test_cpl/k_2/c_1/output/arm_zt.nc','r')
fth=nc.Dataset(clubb_dir+'test_cpl/k_2/c_2/output/arm_zt.nc','r')

# %%
list(range(8,9-1,-1))

# %%
u_r=np.array([1,2,3,4,3,2,1,.5,-.5,-.3,-2,-7])
np.argmin(u_r>0)
usdf=u_r.copy()
usdf[0:4]=0.001
print(usdf[0:4])
u_r[0:4]

# %%
plt.plot(fs['lwp'][:,0,0,0])

# %%
data=fth['thvm'][:,0:150,0,0].T#-ftc['thvm'][:,0:150,0,0].T
#abmax=2#np.max(np.abs(data))
plt.imshow(data,vmin=-abmax,vmax=abmax,cmap='coolwarm',extent=(5,22,0,6))
plt.colorbar()

# %%
plt.plot(ftc['thvm'][460,0:150,0,0],ftc['altitude'][0:150])
plt.plot(fth['thvm'][460,0:150,0,0],ftc['altitude'][0:150])

# %%
days=os.listdir(clubb_dir+'tall2_cpl')
days.sort()
import warnings
warnings.filterwarnings("ignore")

j= 6
for day in days:
    print('.',end='')
    if j == 6:
        #plt.figure(figsize=(15,7))
        fig=plt.figure(figsize=(21,7))
        j=0
    dir2c = clubb_dir+'tall_2c/'+day+'/'
    dir1c = clubb_dir+'tall_1c/'+day+'/'
    dircp = clubb_dir+'tall2_cpl/'+day+'/'
    try:
        fscpc = nc.Dataset(dircp+'k_2/c_1/output/arm_sfc.nc','r')
        fscph = nc.Dataset(dircp+'k_2/c_2/output/arm_sfc.nc','r')
        fc    = nc.Dataset(dircp+'k_2/clusters.nc','r')
        fscp = nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        fs2c = nc.Dataset(dir2c+'k_2/agg_outsfc.nc','r')
        fs1c =nc.Dataset(dir1c+'k_1/c_1/output/arm_sfc.nc','r')
    except Exception as e:
        print(e)
        continue
    #alt=ftcph['altitude'][0:150]

    ax=plt.subplot(2,3,j+1)
    #plt.plot(fs1c['lwp'][:,0,0,0])
    #plt.plot(fs2c['lwp'][:,0,0,0])
    plt.plot(fscpc['lwp'][:,0,0,0])
    plt.plot(fscph['lwp'][:,0,0,0])
    plt.plot(fscp['lwp'][:,0,0,0])
    plt.title(day)
    #data=fc['wbar'][1,:,:150]
    #data[data==0]=float('nan')
    #im=plt.imshow(data.T,origin='lower',cmap='coolwarm',extent=(5,22,0,6))
    #plt.title(day+': '+str(np.max(fc['u_r'][:]))[0:6])
    #cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.01,ax.get_position().height])
    #cbar=plt.colorbar(im,cax=cax)
    #cbar.formatter.set_powerlimits((0, 0))
    
    j=j+1

# %%

# %%
(-3)**3

# %%
import pickle
fp=open('tuneeb.p','rb')
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
    plt.plot(gu,test_line,'b-')
    xaxh=np.nanmax(gu[gu>0])-np.nanmin(gu[gu>0])
    plt.xlim(np.nanmin(gu[gu>0])-xaxh*.1,np.nanmax(gu[gu>0])+xaxh*.3)
    yaxh=np.nanmax(test[test>0])-np.nanmin(test[test>0])
    ymin=np.nanmin(test[test>0])-.1*yaxh
    ymax=np.nanmax(test[test>0])+.3*yaxh
    plt.ylim(ymin,ymax)
    plt.title(name)


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
    plot_test(crs[day],urmx[day],'MAX UR')
    
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

# %%

# %%
day='20170716'
day='20170923'
dirr=clubb_dir+'circ_tune2d/'+day+'/eb_cr0.25k2cc1ccc1/'
dirr=clubb_dir+'tall2_cpl/sgp_20170626/'
ft=nc.Dataset(dirr+'k_2/agg_outzt.nc','r')
fs=nc.Dataset(dirr+'k_2/agg_outsfc.nc','r')
fm=nc.Dataset(dirr+'k_2/agg_outzm.nc','r')
fc=nc.Dataset(dirr+'k_2/clusters.nc','r')
ftc=nc.Dataset(dirr+'k_2/c_1/output/arm_zt.nc','r')
fsc=nc.Dataset(dirr+'k_2/c_1/output/arm_sfc.nc','r')
fmc=nc.Dataset(dirr+'k_2/c_1/output/arm_zm.nc','r')
fth=nc.Dataset(dirr+'k_2/c_2/output/arm_zt.nc','r')
fsh=nc.Dataset(dirr+'k_2/c_2/output/arm_sfc.nc','r')
fmh=nc.Dataset(dirr+'k_2/c_2/output/arm_zm.nc','r')

dirr1=clubb_dir+'circ_tune2d/'+day+'/eb_cr0k1cc0ccc0/'
dirr1=clubb_dir+'tall_1c/sgp_20170626/'
ft1=nc.Dataset(dirr1+'k_1/c_1/output/arm_zt.nc','r')
fs1=nc.Dataset(dirr1+'k_1/c_1/output/arm_sfc.nc','r')
fm1=nc.Dataset(dirr1+'k_1/c_1/output/arm_zm.nc','r')

# %%
ftv=['thlm','p_in_Pa','rtm','thvm','wm','cloud_cover','rtm_forcing','thlm_forcing']
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
    data=ft[v][:,0:150,0,0]
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
    data=fm[v][:,0:150,0,0]
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
t=520
var='wm'
plt.plot(ftc[var][t,0:150,0,0],np.linspace(0,6,150))
plt.plot(fth[var][t,0:150,0,0],np.linspace(0,6,150))

# %%
var='rtm_forcing'
plt.plot(ftc[var][t,0:150,0,0],np.linspace(0,6,150))
plt.plot(fth[var][t,0:150,0,0],np.linspace(0,6,150))

# %%
var='rtm'
plt.plot(ftc[var][t,0:150,0,0],np.linspace(0,6,150))
plt.plot(fth[var][t,0:150,0,0],np.linspace(0,6,150))

# %%
var='thvm'
plt.plot(ftc[var][t,0:150,0,0],np.linspace(0,6,150))
plt.plot(fth[var][t,0:150,0,0],np.linspace(0,6,150))


# %%
# TIMES
# LES is minutes since 2016-01-01 with 10 min dt
# Varanal is seconds since 2012-5-1 with hourly dt
# clubb time is wrong; 0 is 5am local (10am utc) and dt is 1 min
def index_time(date,clb_st=5,les_st=7,var_st=datetime.datetime(2012,5,1,0,0)):
    clb_i=(date.hour-clb_st)*60+date.minute
    if clb_i<0:
        clb_i=clb_i+60*24
    les_i=int((date.hour-les_st)*6+date.minute/10)
    var_i=int((date-var_st).total_seconds()/60/60)+5
    return clb_i,les_i,var_i

# %%
