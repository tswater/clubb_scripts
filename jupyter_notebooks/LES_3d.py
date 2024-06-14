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

# %%

# %%
les_dir='/home/tsw35/tyche/data/LES_FULL/'
les_1c ='/home/tsw35/tyche/data/LES_1C/'
fdir='/home/tsw35/tyche/data/LES_FULL/fr2_20160719_00/'
fdir1chm=les_1c+'trimfr2_20160719_01.nc'
fdir1cht=les_1c+'trimfr2_20160719_00.nc'

# %%
l1chm['time'][1]-l1chm['time'][0]

# %%
l1chm=nc.Dataset(fdir1chm,'r')
l1cht=nc.Dataset(fdir1cht,'r')

# %%
hfiles=os.listdir(fdir)
hfiles.sort()
vpt=np.zeros((16,226,520,520))
H  =np.zeros((16,520,520))
i=0
for file in hfiles:
    print(i)
    lall=nc.Dataset(fdir+file,'r')
    vpt[i,:,:,:]=lall['AVV_THV'][:]
    H[i,:,:]=lall['AVS_SH'][:]
    i=i+1

# %%
l1c['time'][0::6].shape

# %%
clst=np.zeros((520,520))
clst[H[10,:,:]>np.mean(H[10,:,:])]=1

# %%
i=10


# %%
def three_to_one(arr,clst):
    aout=np.zeros((226,2))
    for j in range(226):
        aout[j,0]=np.nanmean(arr[j,:,:][clst==0])
        aout[j,1]=np.nanmean(arr[j,:,:][clst==1])
    return aout


# %%
a=three_to_one(vpt[i,:],clst)

# %%

# %%
vpt_2c=np.zeros((16,226,2))
for i in range(16):
    vpt_2c[i,:]=three_to_one(vpt[i,:],clst)

# %%
diff=l1cht['u2'][1:,0:166]+l1cht['v2'][1:,0:166]-l1chm['u2'][1:,0:166]-l1chm['v2'][1:,0:166]
avg=np.mean(diff,axis=1)
diff=diff/avg[:,None]
vmax=np.percentile(diff,95)
plt.imshow(diff.T,origin='lower',vmin=0,vmax=vmax, cmap='terrain',extent=(7,22,0,10))
plt.xticks(np.arange(7,22,1))
plt.yticks(np.arange(1,10,1),labels=np.arange(.5,5,.5))

# %%
alt=l1cht['kw'][0:166]

# %%
plt.plot(diff[10*6,:],alt)
plt.plot(vpt_2c[10,0:166,1]-320,alt)
plt.plot(vpt_2c[10,0:166,0]-320,alt)

# %%
vptt=vpt[10,:,:,:]

# %%
a=np.argmin(np.abs(vptt-315),axis=0)

# %%
a.shape

# %%
plt.imshow(a.T,origin='lower',cmap='coolwarm')
plt.colorbar()

# %%

# %%
print(a[300,350])
print(a[125,350])

# %%
plt.plot(diff[10*6,:],alt*30)
plt.plot(vpt[10,0:166,300,350]-315,alt*30)
plt.plot(vpt[10,0:166,125,350]-315,alt*30)
plt.xlim(-3.5,2.5)

# %%

# %%
plt.imshow(H[10,:,:].T,origin='lower',cmap='coolwarm')
plt.colorbar()

# %%
days=    [20160625,20160716,20160719,20160720,20170609,
          20170626,20170627,20170629,20170705,20170709,
          20170712,20170716,20170717,20170719,20170720,
          20170728,20170826,20170922,20170923,20170924,
          20180522,20180530,20180618,20180619,20180704,
          20180705,20180523,20180707,20180709,20180710,
          20180711,20180712,20180809,20180811,20180916,
          20180917,20190707,20190709,20190714,20190804,
          20190805]

# %%
vpt_s=np.zeros((len(days),9,166,2))
vpt_a=np.zeros((len(days),9,166,2))
pa_s=np.zeros((len(days),9,166,2))
pa_a=np.zeros((len(days),9,166,2))
circ=np.zeros((len(days),9,166))
j=0
for day in days:
    print(day,end=',')
    fdir='/home/tsw35/tyche/data/LES_FULL/fr2_'+str(day)+'_00/'
    hfiles=os.listdir(fdir)
    hfiles.sort()
    l1cht=nc.Dataset(les_1c+'/trimfr2_'+str(day)+'_01.nc')
    l1chm=nc.Dataset(les_1c+'/trimfr2_'+str(day)+'_00.nc')
    diff=l1cht['u2'][1:,0:166]+l1cht['v2'][1:,0:166]-l1chm['u2'][1:,0:166]-l1chm['v2'][1:,0:166]
    avg=np.mean(diff,axis=1)
    diff=diff/avg[:,None]
    i=0
    for file in hfiles:
        if i<5 or i>13:
            i=i+1
            continue
        lall=nc.Dataset(fdir+file,'r')
        vptt=lall['AVV_THV'][0,:]
        pat =lall['AVV_P'][0,:]
        nm=np.mean(vptt[1,:])+1.5
        a=np.argmin(np.abs(vptt-nm),axis=0)
        minz=np.where(a==np.min(a))
        maxz=np.where(a==np.max(a))
        vpt_a[j,i-5,0:166,1]=np.mean(vptt[0:166,minz[0][:],minz[1][:]],axis=1)
        vpt_a[j,i-5,0:166,0]=np.mean(vptt[0:166,maxz[0][:],maxz[1][:]],axis=1)
        pa_a[j,i-5,0:166,1]=np.mean(pat[0:166,minz[0][:],minz[1][:]],axis=1)
        pa_a[j,i-5,0:166,0]=np.mean(pat[0:166,maxz[0][:],maxz[1][:]],axis=1)
        
        minz=np.where(vptt[1,:]<(np.min(vptt[1,:])+.25))
        maxz=np.where(vptt[1,:]>(np.max(vptt[1,:])-.25))
        vpt_s[j,i-5,0:166,1]=np.mean(vptt[0:166,minz[0][:],minz[1][:]],axis=1)
        vpt_s[j,i-5,0:166,0]=np.mean(vptt[0:166,maxz[0][:],maxz[1][:]],axis=1)
        pa_s[j,i-5,0:166,1]=np.mean(pat[0:166,minz[0][:],minz[1][:]],axis=1)
        pa_s[j,i-5,0:166,0]=np.mean(pat[0:166,maxz[0][:],maxz[1][:]],axis=1)
        circ[j,i-5,0:166]=diff[i*6,0:166]
        print(str(i),end=',')
        i=i+1
    j=j+1
    print()

# %%
fp=open('les_vpt.p','wb')
pickle.dump(vpt_s,fp)
pickle.dump(vpt_a,fp)
pickle.dump(pa_s,fp)
pickle.dump(pa_a,fp)
pickle.dump(circ,fp)
fp.close()

# %%

# %%
dbfile = open('les_vpt.p', 'rb')     
vpt = pickle.load(dbfile)
circ=pickle.load(dbfile)

# %%
vpt.shape

# %%
i=0
j=0
alt=np.linspace(0,5000,166)
plt.figure(figsize=(18,6))
for day in days:
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    plt.subplot(1,6,j+1)
    vpt1=vpt[i,3,:,0]
    vpt2=vpt[i,3,:,1]
    ccc=circ[i,3,:]*2
    mn=min(np.min(vpt1),np.min(vpt2))
    mn2=np.min(vpt1)
    plt.plot(vpt2-mn2,alt)
    plt.plot(vpt1-mn2,alt)
    plt.plot(ccc-np.min(ccc),alt)
    plt.xlim(-1.5,6)
    plt.ylim(0,4000)
    plt.title(day)
    i=i+1
    j=j+1

# %%
i=0
j=0
cdir='/home/tsw35/tyche/clubb/sgp_2c_hf/'
a_l=[]
alt=np.linspace(0,5000,166)
plt.figure(figsize=(18,6))
for day in days:
    fp1=nc.Dataset(cdir+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fp2=nc.Dataset(cdir+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    plt.subplot(1,6,j+1)
    vpt1=vpt[i,0,:,0]
    vpt2=vpt[i,0,:,1]
    alt2=fp1['altitude'][:]
    vptc1=fp1['thvm'][60*(12-5),:,0,0]
    vptc2=fp2['thvm'][60*(12-5),:,0,0]
    #mn=min(np.min(vpt1),np.min(vpt2))
    plt.plot([0,0],[0,5000],'k--',alpha=.2)
    plt.plot(vpt1-vpt2,alt,'r')
    les_alt=alt[np.argmax((vpt1-vpt2)[1:]<0)]
    a_l.append(les_alt)
    plt.scatter([0],les_alt)
    plt.plot(vptc2-vptc1,alt2,'r--')
    plt.scatter([0],alt2[np.argmax((vptc2-vptc1)[1:]<0)])
    plt.xlim(-2,3)
    plt.ylim(0,5000)
    #plt.xlim(min(np.min(vpt2),np.min(vpt1))-1,max(np.max(vpt1),np.max(vpt2))+1)
    plt.title(day)
    i=i+1
    j=j+1

# %%
print(np.mean(a_l))

# %%
i=0
j=0
cdir='/home/tsw35/tyche/clubb/sgp_nocpl_i/'
alt=np.linspace(0,5000,166)
plt.figure(figsize=(18,6))
for day in days:
    fp1=nc.Dataset(cdir+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fp2=nc.Dataset(cdir+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    plt.subplot(1,6,j+1)
    vpt1=vpt[i,5,:,0]
    vpt2=vpt[i,5,:,1]
    alt2=fp1['altitude'][:]
    vptc1=fp1['thvm'][720,:,0,0]
    vptc2=fp2['thvm'][720,:,0,0]
    #mn=min(np.min(vpt1),np.min(vpt2))
    plt.plot(vpt2-vpt2[0],alt,'b')
    plt.plot(vpt1-vpt2[0],alt,'r')
    plt.plot(vptc1-vptc1[0],alt2,'b--')
    plt.plot(vptc2-vptc1[0],alt2,'r--')
    plt.xlim(-1,8)
    plt.ylim(0,5000)
    #plt.xlim(min(np.min(vpt2),np.min(vpt1))-1,max(np.max(vpt1),np.max(vpt2))+1)
    plt.title(day)
    i=i+1
    j=j+1

# %%
np.argmax([])

# %%

# %%
i=0
j=0
cdir='/home/tsw35/tyche/clubb/sgp_nocpl_i/'
alt=np.linspace(0,5000,166)
plt.figure(figsize=(18,6))
for day in days:
    fp1=nc.Dataset(cdir+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fp2=nc.Dataset(cdir+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    plt.subplot(1,6,j+1)
    vpt1=vpt[i,4,:,0]
    vpt2=vpt[i,4,:,1]
    alt2=fp1['altitude'][:]
    vptc1=fp1['thvm'][600,:,0,0]
    vptc2=fp2['thvm'][600,:,0,0]
    #mn=min(np.min(vpt1),np.min(vpt2))
    vptt=(vpt2+vpt1)/2
    vptc=(vptc2+vptc1)/2
    plt.plot(vptt-vptt[0],alt,'r')
    plt.plot(vptc-vptc[0],alt2,'r--')
    plt.ylim(0,5000)
    plt.xlim(-1,20)
    plt.title(day)
    i=i+1
    j=j+1

# %%
fp1['thvm'][:].shape

# %%
a=np.arange(3,21)
a=a.reshape(3,6)
b=np.where(a>12)

# %%
print(a[b[0][:],b[1][:]])

# %%
np.argmin((vptc1-vptc2)<0)


# %%
def min_idx(curv):
    np.argmin(curv<0)
