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
import scipy as sci
import scipy.ndimage
import os
import rasterio
import seaborn as sns
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
sns.set_theme()
plt.rcParams.update({'figure.max_open_warning': 0})

# %%
lwd='/home/tsw35/xTyc_shared/clasp/fr2/fr2_20170626_00/'
lnd='/home/tsw35/xSot_shared/WRFLES/fr2_20170626_00/'
lwdh='/home/tsw35/xTyc_shared/clasp/fr2/fr2_20170626_01/'
bsd='wrfout_d01_2017-06-26_12:00:00'

# %%
flist=os.listdir(lwd)
flist.sort()
u_hmg=np.zeros((13,83,520,520))
u_het=np.zeros((13,83,520,520))
u_nwd=np.zeros((13,83,520,520))

v_hmg=np.zeros((13,83,520,520))
v_het=np.zeros((13,83,520,520))
v_nwd=np.zeros((13,83,520,520))

t_hmg=np.zeros((13,83,520,520))
t_het=np.zeros((13,83,520,520))
t_nwd=np.zeros((13,83,520,520)) #AVV_THV

hfx=np.zeros((13,520,520))

lwp_hmg=np.zeros((13,520,520))
lwp_het=np.zeros((13,520,520))
lwp_nwd=np.zeros((13,520,520)) #AVS_LWP
t=0
for file in flist:
    print(file)
    if '12:00' in file:
        continue
    elif '13:00' in file:
        continue
    elif '14:00' in file:
        continue
    fp1=nc.Dataset(lwdh+file,'r')
    fp2=nc.Dataset(lwd+file,'r')
    fp3=nc.Dataset(lnd+file,'r')
    fp4=nc.Dataset(lnd+'wrfout'+file[4:],'r')
    
    hfx[t,:]=fp2['AVS_SH'][0,:]
    
    u_hmg[t,:]=fp1['AVV_U'][0,0:166:2,:,:]
    u_het[t,:]=fp2['AVV_U'][0,0:166:2,:,:]
    u_nwd[t,:]=fp4['U'][0,0:166:2,:,:-1]
    
    v_hmg[t,:]=fp1['AVV_V'][0,0:166:2,:,:]
    v_het[t,:]=fp2['AVV_V'][0,0:166:2,:,:]
    v_nwd[t,:]=fp4['V'][0,0:166:2,:-1,:]
    
    t_hmg[t,:]=fp1['AVV_THV'][0,0:166:2,:,:]
    t_het[t,:]=fp2['AVV_THV'][0,0:166:2,:,:]
    t_nwd[t,:]=fp3['AVV_THV'][0,0:166:2,:,:]
    
    t_hmg[t,t_hmg[t,:]<np.percentile(t_hmg[t,:],.1)]=np.percentile(t_hmg[t,:],.1)
    t_het[t,t_hmg[t,:]<np.percentile(t_het[t,:],.1)]=np.percentile(t_het[t,:],.1)
    t_nwd[t,t_hmg[t,:]<np.percentile(t_nwd[t,:],.1)]=np.percentile(t_nwd[t,:],.1)
    
    lwp_hmg[t,:]=fp1['AVS_LWP'][0,:]
    lwp_het[t,:]=fp2['AVS_LWP'][0,:]
    lwp_nwd[t,:]=fp3['AVS_LWP'][0,:]

    t=t+1

# %%
data1=sci.ndimage.filters.gaussian_filter(v_het[4,20,:],15,mode='reflect')
data2=sci.ndimage.filters.gaussian_filter(v_nwd[4,20,:],15,mode='reflect')
vmax=max(np.max(np.abs(data1)),np.max(np.abs(data2)))
vmin=-vmax
plt.subplot(1,2,1)
plt.imshow(data1,origin='lower',cmap='PRGn',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.title('Background Wind')
plt.axis(False)
plt.subplot(1,2,2)
plt.imshow(data2,origin='lower',cmap='PRGn',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.axis(False)
plt.title('No Background Wind')
plt.show()

# %%
plt.plot(v_het[:,10,200,200])
plt.show()

# %%
data1=sci.ndimage.filters.gaussian_filter(v_het[4,10,:],15,mode='reflect')
data1=data1-np.mean(data1)
data2=sci.ndimage.filters.gaussian_filter(v_nwd[4,10,:],15,mode='reflect')
data2=data2-np.mean(data2)
vmax=max(np.max(np.abs(data1)),np.max(np.abs(data2)))
vmin=-vmax
plt.subplot(1,2,1)
plt.imshow(data1,origin='lower',cmap='PRGn',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.title('Background Wind')
plt.axis(False)
plt.subplot(1,2,2)
plt.imshow(data2,origin='lower',cmap='PRGn',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.axis(False)
plt.title('No Background Wind')
plt.show()

# %%
div=np.gradient(sci.ndimage.filters.gaussian_filter(u_het[4,10,:],15,mode='reflect'))[1]/3+np.gradient(sci.ndimage.filters.gaussian_filter(v_het[4,10,:],15,mode='reflect'))[0]/3

# %%
div_nwd=np.gradient(sci.ndimage.filters.gaussian_filter(u_nwd[4,10,:],15,mode='reflect'))[1]/3+np.gradient(sci.ndimage.filters.gaussian_filter(v_nwd[4,10,:],15,mode='reflect'))[0]/3

# %%
data1=div
data2=div_nwd
vmax=max(np.max(np.abs(data1)),np.max(np.abs(data2)))
vmin=-vmax
plt.subplot(1,2,1)
plt.imshow(data1,origin='lower',cmap='PRGn',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.title('Background Wind')
plt.axis(False)
plt.subplot(1,2,2)
plt.imshow(data2,origin='lower',cmap='PRGn',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.axis(False)
plt.title('No Background Wind')
plt.show()

# %%
data1=t_het[4,10,:]
data2=t_nwd[4,10,:]
vmax=306
vmin=304
plt.subplot(1,2,1)
plt.imshow(data1,origin='lower',cmap='coolwarm',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.title('Wind '+str(np.std(data1))[0:7])
plt.axis(False)
plt.subplot(1,2,2)
plt.imshow(data2,origin='lower',cmap='coolwarm',vmin=vmin,vmax=vmax)
plt.colorbar(shrink=0.5)
plt.axis(False)
plt.title('No Wind '+str(np.std(data2))[0:7])
plt.show()
print(np.std(data1)/np.std(data2))

# %%
fig=plt.figure()
fig.patch.set_facecolor('lightgrey')
data1=lwp_het[4,:,:]
data2=lwp_nwd[4,:,:]
data1[data1==0]=float('nan')
data2[data2==0]=float('nan')
plt.subplot(1,2,1)
plt.imshow(data1,origin='lower',cmap='Blues',vmin=0,vmax=.25)
plt.colorbar(shrink=0.5)
plt.title('Wind '+str(np.std(data1))[0:7])
plt.axis(False)
plt.subplot(1,2,2)
plt.imshow(data2,origin='lower',cmap='Blues',vmin=0,vmax=.25)
plt.colorbar(shrink=0.5)
plt.axis(False)
plt.title('No Wind '+str(np.std(data2))[0:7])
plt.show()
print(np.std(data1)/np.std(data2))

# %%
np.std(data1)/np.std(data2)

# %%
t_het_std=np.zeros((13,83))
t_nwd_std=np.zeros((13,83))
d_het=np.zeros((13,83))
d_nwd=np.zeros((13,83))
for t in range(13):
    print('.',end='',flush=True)
    for i in range(83):
        t_het_std[t,i]=np.std(t_het[t,i,40:-40,40:-40])
        t_nwd_std[t,i]=np.std(t_nwd[t,i,40:-40,40:-40])
        d_het[t,i]=np.mean(np.abs(np.gradient(sci.ndimage.filters.gaussian_filter(u_het[t,i,:],10,mode='reflect'))[1]/3+np.gradient(sci.ndimage.filters.gaussian_filter(v_het[t,i,:],10,mode='reflect'))[0]/3))
        d_nwd[t,i]=np.mean(np.abs(np.gradient(sci.ndimage.filters.gaussian_filter(u_nwd[t,i,:],10,mode='reflect'))[1]/3+np.gradient(sci.ndimage.filters.gaussian_filter(v_nwd[t,i,:],10,mode='reflect'))[0]/3))

# %%
for t in range(4,12):
    plt.figure()
    plt.plot(d_het[t,:]/t_het_std[t,:],np.linspace(0,2500,83))
    plt.plot(d_nwd[t,:]/t_nwd_std[t,:],np.linspace(0,2500,83))
    plt.ylim(0,1500)
    plt.xlim(0,.1)
    plt.legend([r"$\Nabla$ $\sigma_T$"])
plt.show()

# %%
plt.imshow()

# %%
print(np.mean(np.abs(div)))

# %%
np.mean(np.abs(div))/np.mean(np.abs(div_nwd))

# %%
np.std(data1[data1>300])/np.std(data2[data2>300])

# %%
print(np.mean(np.abs(div_nwd)))

# %%
print(np.std(data1[data1>300]))

# %%
print(np.std(data2[data2>300]))

# %%
plt.imshow(data1-data2)
plt.colorbar()
plt.show()

# %%
fp=nc.Dataset('/home/tsw35/xSot_shared/WRFLES/fr2_20170626_00/wrfout_d01_2017-06-26_22:00:00')

# %%
data=sci.ndimage.filters.gaussian_filter(fp['V'][0,75,:,:],15,mode='reflect')
#data=u_het[3,5,:,:]
plt.imshow(data,origin='lower',cmap='PRGn')
plt.colorbar()
plt.show()

# %%
flist=os.listdir(lwd)
flist.sort()
#flist[3]
fp2=nc.Dataset(lwd+flist[10],'r')
data=sci.ndimage.filters.gaussian_filter(fp2['AVV_V'][0,75,:,:],15,mode='reflect')
#data=u_het[3,5,:,:]
plt.imshow(data,origin='lower',cmap='PRGn')
plt.colorbar()
plt.show()

# %%
u90=np.zeros((166,))
u10=np.zeros((166,))
v90=np.zeros((166,))
v10=np.zeros((166,))

u90n=np.zeros((166,))
u10n=np.zeros((166,))
v90n=np.zeros((166,))
v10n=np.zeros((166,))

dt=np.zeros((166,))
dtn=np.zeros((166,))

for i in range(166):
    print('.',end='')
    u=sci.ndimage.filters.gaussian_filter(fp2['AVV_U'][0,i,:,:],15,mode='reflect')
    v=sci.ndimage.filters.gaussian_filter(fp2['AVV_V'][0,i,:,:],15,mode='reflect')
    un=sci.ndimage.filters.gaussian_filter(fp['U'][0,i,:,:],15,mode='reflect')
    vn=sci.ndimage.filters.gaussian_filter(fp['V'][0,i,:,:],15,mode='reflect')
    t=sci.ndimage.filters.gaussian_filter(fp2['AVV_TH'][0,i,:,:],10,mode='reflect')
    tn=sci.ndimage.filters.gaussian_filter(fp['T'][0,i,:,:],10,mode='reflect')

    u90[i]=np.percentile(u,90)
    u10[i]=np.percentile(u,10)
    v90[i]=np.percentile(v,90)
    v10[i]=np.percentile(v,10)

    u90n[i]=np.percentile(un,90)
    u10n[i]=np.percentile(un,10)
    v90n[i]=np.percentile(vn,90)
    v10n[i]=np.percentile(vn,10)
    
    t1=t[t>np.percentile(t,75)]
    t2=t[t<np.percentile(t,25)]
    
    t1n=tn[tn>np.percentile(tn,75)]
    t2n=tn[tn<np.percentile(tn,25)]
    
    dt[i]=np.mean(t1)-np.mean(t2)
    dtn[i]=np.mean(t1n)-np.mean(t2n)

# %%
uu=(u90-u10)/2
vv=(v90-v10)/2
uun=(u90n-u10n)/2
vvn=(v90n-v10n)/2
u2=np.sqrt(uu**2+vv**2)
u2n=np.sqrt(uun**2+vvn**2)
 
um=np.mean(fp2['AVV_U'][0,0:166,:,:],axis=(1,2))
vm=np.mean(fp2['AVV_V'][0,0:166,:,:],axis=(1,2))
umn=np.mean(fp['U'][0,0:166,:,:],axis=(1,2))
vmn=np.mean(fp['V'][0,0:166,:,:],axis=(1,2))

# %%

# %%
plt.plot(u2,np.linspace(0,165,166))
plt.plot(u2n,np.linspace(0,165,166))
plt.show()

# %%
plt.plot(dt,np.linspace(0,165,166))
plt.plot(dtn,np.linspace(0,165,166))
plt.show()

# %%
plt.plot(np.sqrt(um**2+vm**2),np.linspace(0,165,166))
plt.plot(np.sqrt(umn**2+vmn**2),np.linspace(0,165,166))
plt.show()

# %%
plt.plot(u2n/u2,np.linspace(0,5,166))
plt.plot(dtn/dt,np.linspace(0,5,166))
#plt.plot(u2n/u2-dtn/dt,np.linspace(0,5,166))
#plt.plot(np.sqrt(vm**2)-1,np.linspace(0,5,166))
plt.ylim(0,3)
plt.xlim(0,3)
plt.show()

# %%
plt.plot(dt[0:50],np.sqrt(u90**2+v10**2)[0:50])
plt.plot(dtn[0:50],np.sqrt(u90n**2+v10n**2)[0:50])
plt.show()

# %%
np.mean(dtn[0:50]/dt[0:50])

# %%
np.mean(u2n[0:50]/u2[0:50])

# %%
