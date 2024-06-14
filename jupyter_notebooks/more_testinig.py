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
import os
import seaborn as sns
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
sns.set_theme()
plt.rcParams.update({'figure.max_open_warning': 0})

# %% language="bash"
# ls /home/tsw35/tyche/clubb/test/test1.00_1_1_0/k_1/c_1/output

# %%
tdir='/home/tsw35/tyche/clubb/test/'
clubb_dir='/home/tsw35/tyche/clubb/'
dir1c=tdir+'test1.00_1_1_0'
dir2c=tdir+'test1.00_1_2_0'
#fps1c=nc.Dataset(dir1c+'/k_1/c_1/output/arm_sfc.nc','r')
fps2c=nc.Dataset(dir2c+'/k_2/agg_outsfc.nc','r')
fpz2clo=nc.Dataset(dir2c+'/k_2/c_1/output/arm_zt.nc','r')
fpz2chi=nc.Dataset(dir2c+'/k_2/c_2/output/arm_zt.nc','r')

# %%
plt.figure(figsize=(8,4))
t=np.linspace(5,22,fps2c['lwp'][:,0,0,0].shape[0])
#plt.plot(t,fps1c['lwp'][:,0,0,0],'k-',linewidth=1)
plt.plot(t,fps2c['lwp'][:,0,0,0],'k--')
leglist=['2C']
for folder in os.listdir(tdir):
    if '_0' in folder:
        continue
    #if '10' not in folder:
     #   continue
    leglist.append(folder[4:11])
    fps=nc.Dataset(tdir+folder+'/k_2/agg_outsfc.nc','r')
    lwp=fps['lwp'][0:,0,0,0]
    plt.plot(t[0:lwp.shape[0]],lwp)
plt.legend(leglist)

# %%
plt.figure(figsize=(8,4))
t=np.linspace(5,22,int(fps2c['lwp'][:,0,0,0].shape[0]/5)+1)
#plt.plot(t,fps1c['lwp'][:,0,0,0],'k-',linewidth=1)
#plt.plot(t,fps2c['lwp'][:,0,0,0],'k--')
leglist=[]
for folder in os.listdir(tdir):
    if '_0' in folder:
        continue
    try:
        fps=nc.Dataset(tdir+folder+'/k_2/clusters.nc','r')
        lwp=fps['z_circ'][:]
    except:
        continue
    #urz=fps['ur_z'][:,0]
    leglist.append(folder[4:11])
    plt.plot(t[0:lwp.shape[0]],lwp)
plt.legend(leglist)

# %%
plt.plot(lwp)

# %%
import warnings
warnings.filterwarnings('ignore')
plt.figure(figsize=(16,8))
alt=fpz2clo['altitude'][:]
t=590
var='rtm'
lo2c=fpz2clo[var][t,:,0,0]
hi2c=fpz2chi[var][t,:,0,0]
d2c=hi2c-lo2c
leglist=[]

for folder in os.listdir(tdir):
    if '_0' in folder:
        continue
    print(tdir+folder,end='   ')
    fpzlo=nc.Dataset(tdir+folder+'/k_2/c_1/output/arm_zt.nc','r')
    fpzhi=nc.Dataset(tdir+folder+'/k_2/c_2/output/arm_zt.nc','r')
    try:
        plt.subplot(141)
        plt.plot(fpzlo[var][t,:,0,0],alt)
        plt.ylim(500,4500)
        #plt.xlim(315,330)
        plt.subplot(142)
        plt.plot(fpzhi[var][t,:,0,0],alt)
        plt.ylim(500,4500)
        #plt.xlim(315,330)
        plt.subplot(143)
        plt.plot(fpzhi[var][t,:,0,0]-fpzlo[var][t,:,0,0],alt)
        plt.ylim(500,4500)
        plt.subplot(144)
        plt.plot([0,0],[0,0])
    except:
        continue
    leglist.append(folder[4:11])
plt.subplot(144)
plt.legend(leglist)
    

# %%
#### FORCINGS ####

# %%
fc10_l = read_forcings('/home/tsw35/tyche/clubb/test/test1.00_10_2_1/k_2/c_1/arm_forcings_f.in')
fc1_l = read_forcings('/home/tsw35/tyche/clubb/test/test1.00_1_2_1/k_2/c_1/arm_forcings_f.in')
fc0_l = read_forcings('/home/tsw35/tyche/clubb/test/test1.00_1_2_0/k_2/c_1/arm_forcings_f.in')
fc10_h = read_forcings('/home/tsw35/tyche/clubb/test/test1.00_10_2_1/k_2/c_2/arm_forcings_f.in')
fc1_h = read_forcings('/home/tsw35/tyche/clubb/test/test1.00_1_2_1/k_2/c_2/arm_forcings_f.in')
fc0_h = read_forcings('/home/tsw35/tyche/clubb/test/test1.00_1_2_0/k_2/c_2/arm_forcings_f.in')

# %%
abmax=max(np.max(np.abs(data1)),np.max(np.abs(data2)))
plt.subplot(121)
plt.imshow(data1,vmin=-abmax,vmax=abmax,cmap='coolwarm',origin='lower',extent=[5,22,0,12])

# %%
var='rtm_f[kg\kg\s]' #'T_f[K\s]' #'rtm_f[kg\kg\s]'
plt.figure(figsize=(15,5))

#### 10 case ####
data1=fc10_l[var].T-fc0_l[var].T
data2=fc10_h[var].T-fc0_h[var].T
abmax=max(np.max(np.abs(data1)),np.max(np.abs(data2)))
plt.subplot(2,2,1)
plt.imshow(data1,vmin=-abmax,vmax=abmax,cmap='coolwarm',origin='lower',extent=[5,22,0,12])
plt.title('LOW -- 10w')
plt.ylim(0,2)
plt.xlim(10,20)
plt.subplot(2,2,2)
plt.imshow(data2,vmin=-abmax,vmax=abmax,cmap='coolwarm',origin='lower',extent=[5,22,0,12])
plt.title('HIGH -- 10w')
plt.ylim(0,2)
plt.xlim(10,20)

#### 1 case ####
data1=fc10_l[var].T-fc0_l[var].T
data2=fc10_h[var].T-fc0_h[var].T
abmax=max(np.max(np.abs(data1)),np.max(np.abs(data2)))
plt.subplot(2,2,3)
plt.imshow(data1,vmin=-abmax,vmax=abmax,cmap='coolwarm',origin='lower',extent=[5,22,0,12])
plt.title('LOW -- 1w')
plt.ylim(0,2)
plt.xlim(10,20)
plt.subplot(2,2,4)
plt.imshow(data2,vmin=-abmax,vmax=abmax,cmap='coolwarm',origin='lower',extent=[5,22,0,12])
plt.title('HIGH -- 1w')
plt.ylim(0,2)
plt.xlim(10,20)

# %%
fpzthi=nc.Dataset('/home/tsw35/tyche/clubb/test/test1.00_1_2_1/k_2/c_2/output/arm_zt.nc','r')
fpztlo=nc.Dataset('/home/tsw35/tyche/clubb/test/test1.00_1_2_1/k_2/c_1/output/arm_zt.nc','r')
fpzmhi=nc.Dataset('/home/tsw35/tyche/clubb/test/test1.00_1_2_1/k_2/c_2/output/arm_zm.nc','r')
fpzmlo=nc.Dataset('/home/tsw35/tyche/clubb/test/test1.00_1_2_1/k_2/c_1/output/arm_zm.nc','r')
alt=fpzthi['altitude'][:]/1000
alt.shape

# %%

# %%

# %%
fig,[ax1,ax2,ax3,ax4]=plt.subplots(1,4,figsize=(12,10))
times=np.linspace(0,1000,100)
t_title=np.linspace(5,22,100)
alt=fpzthi['altitude'][0:175]/1000

data_hi=fpzthi['thvm'][times.astype(int),0:175,0,0]#-mean
data_lo=fpztlo['thvm'][times.astype(int),0:175,0,0]#-mean
T_hi=fpzthi['T_in_K'][times.astype(int),0:175,0,0]
T_lo=fpztlo['T_in_K'][times.astype(int),0:175,0,0]
rtm_hi=fpzthi['rtm'][times.astype(int),0:175,0,0]
rtm_lo=fpztlo['rtm'][times.astype(int),0:175,0,0]
ri_hi=fpzmhi['Richardson_num'][times.astype(int),0:175,0,0]
ri_lo=fpzmlo['Richardson_num'][times.astype(int),0:175,0,0]

def animate(i):
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    
    ax1.plot(data_hi[i,:],alt,c='r')
    ax1.plot(data_lo[i,:],alt,c='b')
    ax1.set_title("VPT")
    ax1.set_xlim(300,330)
    
    ax2.plot(T_hi[i,:],alt,c='r')
    ax2.plot(T_lo[i,:],alt,c='b')
    ax2.set_title("T_in_K")
    ax2.set_xlim(260,310)
    
    ax3.plot(rtm_hi[i,:],alt,c='r')
    ax3.plot(rtm_lo[i,:],alt,c='b')
    ax3.set_title("rtm")
    ax3.set_xlim(-.002,.0175)
    
    ax4.plot(ri_hi[i,:],alt,c='r')
    ax4.plot(ri_lo[i,:],alt,c='b')
    ax4.set_title("Richardson Number")
    ax4.set_xlim(-1,20)
    
    plt.suptitle(t_title[i])

    return fig
    

ani=FuncAnimation(fig,animate,frames=100,interval=100,repeat=True)
#FFwriter = animation.FFMpegWriter(fps=10)
#ani.save('ani_1.mp4',writer=FFwriter)
HTML(ani.to_jshtml())

# %%
grad=np.gradient(data_hi.T,axis=1)
plt.figure(figsize=(12,8))
plt.imshow(data_hi-np.mean(data_hi,axis=0),origin='lower',extent=(5,22,0,7))
plt.colorbar()

# %%
fps2c['lwp'].shape


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
fp=nc.Dataset('/home/tsw35/soteria/clubb/data/sgp60varanarap_2012-2019.nc','r')

# %%
for var in fp.variables:
    print(var)

# %%
plt.hist(fp['p_srf_aver'][:],bins=50)

# %%
plt.plot(fp['p_srf_aver'][:])

# %%
plt.plot(fp['lev'][:])

# %%
for var in fp.variables:
    print(var)

# %%
clubb_dir = '/home/tsw35/tyche/clubb/'
fdir = clubb_dir+'sgp_1c'
les_dir= '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
day_list2=os.listdir(fdir)
i=1

fp=nc.Dataset(les_dir+'trimdoz_'+str(day[4:12])+'_01.nc','r')
alt=fp['k'][:]*30
dt1 = datetime.datetime(2016,1,1,0,0)
times_js = []
for t in range(len(fp['time'][:])):
    times_js.append(dt1+datetime.timedelta(minutes=int(fp['time'][t]))-datetime.timedelta(hours=5))
    
plt.figure(figsize=(16,10))
for day in day_list2:
    asdf
    if i==17:
        i=1
        plt.figure(figsize=(16,10))
    plt.subplot(4,4,i)
    i = i+1
    fp=nc.Dataset(les_dir+'trimdoz_'+str(day[4:12])+'_01.nc','r')
    thv=fp['thv'][75,:]
    plt.plot(thv,alt)
    plt.xlim(np.min(thv)-1,np.max(thv[0:100]))
    plt.ylim(0,3000)
    

# %%
clubb_dir = '/home/tsw35/tyche/clubb/'
fdir = clubb_dir+'sgp_1c'
les_dir= '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
day_list2=os.listdir(fdir)
i=1
plt.figure(figsize=(16,10))

for day in day_list2:
    if i==17:
        i=1
        plt.figure(figsize=(16,10))
    plt.subplot(4,4,i)
    i = i+1
    try:
        fpc1=nc.Dataset(clubb_dir+'sgp_nocpl_m/'+day+'/k_2/c_1/output/arm_zt.nc','r')
        fpc2=nc.Dataset(clubb_dir+'sgp_nocpl_m/'+day+'/k_2/c_2/output/arm_zt.nc','r')
    except Exception as e:
        print(day)
        print(e)
        print()
    alt=fpc1['altitude'][:]
    thv1=fpc1['thlm'][700,:,0,0]
    thv2=fpc2['thlm'][700,:,0,0]
    plt.plot(thv1,alt)
    plt.plot(thv2,alt)
    #plt.xlim(min(np.min(thv1),np.min(thv2))-1,max(np.min(thv1),np.min(thv2))+5)
    plt.ylim(0,6000)

# %%
print(day)

# %%
int_days=[20160625,20160716,20160719,20160720,20170609,
          20170626,20170627,20170629,20170705,20170709,
          20170712,20170716,20170717,20170719,20170720,
          20170728,20170826,20170922,20170923,20170924,
          20180522,20180530,20180618,20180619,20180704,
          20180705,20180523,20180707,20180709,20180710,
          20180711,20180712,20180809,20180811,20180916,
          20180917,20190707,20190709,20190714,20190804,
          20190805]
int_days.sort()
days_nm = []
for i in int_days:
    e='sgp_'+str(i)
    days_nm.append(e)

# %%

# %%
clubb_dir = '/home/tsw35/tyche/clubb/'
fdir = clubb_dir+'sgp_1c'
les_dir= '/home/tsw35/tyche/data/LES_1C/'
day=days_nm[4]
i=1

fp=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
alt=fp['k'][:]*30
dt1 = datetime.datetime(2016,1,1,0,0)
times_js = []
for t in range(len(fp['time'][:])):
    times_js.append(dt1+datetime.timedelta(minutes=int(fp['time'][t]))-datetime.timedelta(hours=5))
    
plt.figure(figsize=(24,10))
for day in days_nm[:]:
    if i==17:
        #break
        i=1
        plt.figure(figsize=(24,10))
    plt.subplot(4,4,i)
    i = i+1
    fphm=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_01.nc','r')
    fpht=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
    #plt.plot(fp['tke'][:]-fp2['tke'][:])
    #diff=np.gradient(fpht['thv'][1:,0:134],axis=1)#-fphm['v2'][1:,0:134]
    diff=fpht['qv'][1:,0:134]-fphm['qv'][1:,0:134]
    avg=np.mean(diff,axis=1)
    #diff=diff/avg[:,None]
    #diff=np.gradient(diff,axis=1)
    #abmax=np.max(np.abs(diff))
    abmax=np.percentile(np.abs(diff),98)
    plt.imshow(diff.T,origin='lower',vmin=-abmax,vmax=abmax,cmap='coolwarm',extent=(7,22,0,4))
    #plt.imshow(diff.T,origin='lower',vmin=0,vmax=abmax,cmap='terrain',extent=(7,22,0,4))
    #plt.imshow(diff.T,origin='lower',vmin=-.1,vmax=.3,cmap='terrain',extent=(7,22,0,4))
    plt.colorbar()
    plt.title(day)
    #plt.xlim(np.min(thv),np.max(var[0:135]))
    #plt.ylim(0,4000)

# %%
les_dir= '/home/tsw35/tyche/data/LES_1C/'
fp=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
alt=fp['k'][:]*30
dt1 = datetime.datetime(2016,1,1,0,0)
times_js = []
for t in range(len(fp['time'][:])):
    times_js.append(dt1+datetime.timedelta(minutes=int(fp['time'][t]))-datetime.timedelta(hours=5))
days_nm.sort()
for day in days_nm[:]:
    fphm=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_01.nc','r')
    fpht=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
    plt.figure(figsize=(18,8))
    
    plt.subplot(1,2,1)
    #diff=fpht['v'][1:,0:166]+fpht['u'][1:,0:166]-fphm['v'][1:,0:166]-fphm['u'][1:,0:166]
    diff=fpht['v2'][1:,0:166]+fpht['u2'][1:,0:166]-fphm['v2'][1:,0:166]-fphm['u2'][1:,0:166]
    #diff = fpht['w'][1:,0:166]-fphm['w'][1:,0:166]
    stren=np.mean(fpht['wth'][:,0]) # strength of circulation
    plt.title(day+'   str: %.3f' %stren)
    avg=np.mean(diff,axis=1)
    diff=diff/avg[:,None]
    #vmin=-abmax,vmax=abmax
    vmax=np.percentile(np.abs(diff),95)
    plt.imshow(diff.T,origin='lower',vmin=0,vmax=vmax, cmap='terrain',extent=(7,22,0,10))
    #plt.imshow(diff.T,origin='lower',vmin=-vmax,vmax=vmax, cmap='terrain',extent=(7,22,0,10))
    plt.xticks(np.arange(7,22,1))
    plt.yticks(np.arange(1,10,1),labels=np.arange(.5,5,.5))
    plt.colorbar()
    
    thv=np.gradient(fpht['thv'][:,0:166],axis=1)
    
    
    plt.subplot(1,2,2)
    diff = fpht['w2'][1:,0:166]-fphm['w2'][1:,0:166]
    vmax=np.percentile(np.abs(diff),95)
    plt.imshow(diff.T,origin='lower',vmin=-vmax,vmax=vmax, cmap='terrain',extent=(7,22,0,10))
    plt.colorbar()
    plt.title(day+'   str: %.3f' %stren)
    #plt.scatter(xx,yy,c='k')
    plt.xticks(np.arange(7,22,1))
    plt.yticks(np.arange(1,10,1),labels=np.arange(.5,5,.5))

# %%
import pickle
fp=open('les_vpt.p','rb')
vpt_s=pickle.load(fp)
vpt_a=pickle.load(fp)
pa_s=pickle.load(fp)
pa_a=pickle.load(fp)
circ=pickle.load(fp)
fp.close()

# %%

# %%
les_dir= '/home/tsw35/tyche/data/LES_1C/'
fp=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
alt=fp['k'][:]*30
dt1 = datetime.datetime(2016,1,1,0,0)
times_js = []
for t in range(len(fp['time'][:])):
    times_js.append(dt1+datetime.timedelta(minutes=int(fp['time'][t]))-datetime.timedelta(hours=5))

i=0
for day in days_nm:
    fphm=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_01.nc','r')
    fpht=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
    plt.figure()
    
    diff=fpht['v2'][1:,0:166]+fpht['u2'][1:,0:166]-fphm['v2'][1:,0:166]-fphm['u2'][1:,0:166]
    #diff=fpht['qv2'][1:,0:166]-fphm['qv2'][1:,0:166]
    print(np.mean(diff[:,33:166]))
    stren=np.mean(fpht['wth'][:,0]) # strength of circulation
    plt.title(day+'   str: %.3f' %stren)
    avg=np.mean(diff,axis=1)
    diff=diff/avg[:,None]
    vmax=np.percentile(np.abs(diff),95)
    plt.imshow(diff.T,origin='lower',vmin=0,vmax=vmax, cmap='terrain',extent=(7,22,0,10))
    plt.xticks(np.arange(7,22,1))
    plt.yticks(np.arange(1,10,1),labels=np.arange(.5,5,.5))
    plt.colorbar(label=r"$\sigma_u^2$ (Circulation)")
    
    diff=vpt_a[i,:,:,0]-vpt_a[i,:,:,1]
    #a,b=getlinecrs(diff)
    a,b,c=getlinehdif(vpt_a[i,:,:,0],vpt_a[i,:,:,1],2.5)
    #a,b,c=getlinepa(vpt_a[i,:,:,0],vpt_a[i,:,:,1],pa_a[i,:,:,0],pa_a[i,:,:,1],2)
    plt.plot([12,13,14,15,16,17,18,19,20],a,'k-.')
    plt.plot([12,13,14,15,16,17,18,19,20],b,'k--*')
    plt.plot([12,13,14,15,16,17,18,19,20],c,'k--o')
    plt.legend([r'$z_{hot}$ at $\theta_{v max}$',r'$z_{cold}$ at $\theta_{v max}$',r'$z$ at $\theta_{v crit}$'],framealpha=0.99)
    plt.xlabel('Time (hr)')
    plt.ylabel('Elevation')
    i=i+1

# %%
les_dir= '/home/tsw35/tyche/data/LES_1C/'
fp=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
alt=fp['k'][:]*30
dt1 = datetime.datetime(2016,1,1,0,0)
times_js = []
for t in range(len(fp['time'][:])):
    times_js.append(dt1+datetime.timedelta(minutes=int(fp['time'][t]))-datetime.timedelta(hours=5))

i=0
for day in ['sgp_20160720']:
    fphm=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_01.nc','r')
    fpht=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
    plt.figure(figsize=(9,4),dpi=300)
    
    diff=fpht['v2'][1:,0:166]+fpht['u2'][1:,0:166]-fphm['v2'][1:,0:166]-fphm['u2'][1:,0:166]
    #diff=fpht['qv2'][1:,0:166]-fphm['qv2'][1:,0:166]
    print(np.mean(diff[:,33:166]))
    stren=np.mean(fpht['wth'][:,0]) # strength of circulation
    #plt.title(day+'   str: %.3f' %stren)
    avg=np.mean(diff,axis=1)
    diff=diff/avg[:,None]
    vmax=np.percentile(np.abs(diff),95)
    im = plt.imshow(diff.T,origin='lower',vmin=0,vmax=vmax, cmap='terrain',extent=(7,22,0,10))
    plt.xticks(np.arange(7,26,4),fontsize=18)
    plt.yticks(np.arange(0,15,5),labels=np.arange(0,7.5,2.5),fontsize=18)
    cb=plt.colorbar()
    cb.set_label(label=r"$\Delta$ $\sigma_u^2/u_0$ (Circulation)",size=18,weight='bold')
    cb.ax.tick_params(labelsize=18)
    i=5
    diff=vpt_a[i,:,:,0]-vpt_a[i,:,:,1]
    #a,b=getlinecrs(diff)
    a,b,c=getlinehdif(vpt_a[i,:,:,0],vpt_a[i,:,:,1],2)
    #a,b,c=getlinepa(vpt_a[i,:,:,0],vpt_a[i,:,:,1],pa_a[i,:,:,0],pa_a[i,:,:,1],2)
    plt.plot([12,13,14,15,16,17,18,19,20],a,'k-.')
    plt.plot([12,13,14,15,16,17,18,19,20],b,'k--*')
    plt.plot([12,13,14,15,16,17,18,19,20],c,'k--o')
    plt.legend([r'$z_{hot}$ of $\theta_{v max}$',r'$z_{cold}$ of $\theta_{v max}$',r'$z$ of $\theta_{v crit}$'],framealpha=0.99,fontsize=18)
    plt.xlabel('Time (hr)',fontsize=20,weight='bold')
    plt.ylabel('Elevation (km)',fontsize=20,weight='bold')
    i=i+1

# %%
les_dir= '/home/tsw35/tyche/data/LES_1C/'
fp=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
alt=fp['k'][:]*30
dt1 = datetime.datetime(2016,1,1,0,0)
times_js = []
for t in range(len(fp['time'][:])):
    times_js.append(dt1+datetime.timedelta(minutes=int(fp['time'][t]))-datetime.timedelta(hours=5))

i=0
for day in days_nm[:]:
    fphm=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_01.nc','r')
    fpht=nc.Dataset(les_dir+'trimfr2_'+str(day[4:12])+'_00.nc','r')
    plt.figure()
    
    #diff=fpht['v'][1:,0:166]+fpht['u'][1:,0:166]-fphm['v'][1:,0:166]-fphm['u'][1:,0:166]
    
    diff=np.sqrt(fpht['v'][1:,0:166]**2+fpht['u'][1:,0:166]**2)-np.sqrt(fphm['v'][1:,0:166]**2+fphm['u'][1:,0:166]**2)
    #diff=fpht['qv2'][1:,0:166]-fphm['qv2'][1:,0:166]
    stren=np.mean(fpht['wth'][:,0]) # strength of circulation
    plt.title(day+'   str: %.3f' %stren)
    #avg=np.mean(diff,axis=1)
    #diff=diff/avg[:,None]
    vmax=np.percentile(np.abs(diff),95)
    plt.imshow(np.abs(diff.T),origin='lower',vmin=0,vmax=vmax, cmap='terrain',extent=(7,22,0,10))
    plt.xticks(np.arange(7,22,1))
    plt.yticks(np.arange(1,10,1),labels=np.arange(.5,5,.5))
    plt.colorbar()
    
    diff=vpt_a[i,:,:,0]-vpt_a[i,:,:,1]
    #a,b=getlinecrs(diff)
    a,b,c=getlinehdif(vpt_a[i,:,:,0],vpt_a[i,:,:,1],2)
    #a,b,c=getlinepa(vpt_a[i,:,:,0],vpt_a[i,:,:,1],pa_a[i,:,:,0],pa_a[i,:,:,1],2)
    plt.plot([12,13,14,15,16,17,18,19,20],a,'k-.')
    plt.plot([12,13,14,15,16,17,18,19,20],b,'k--*')
    plt.plot([12,13,14,15,16,17,18,19,20],c,'k--o')
    
    i=i+1


# %%
def getlinecrs(diff_):
    diff_[:,0:20]=0
    out1=np.zeros((9,))
    out2=np.zeros((9,))
    for j in range(9):
        try:
            out1[j]=np.argmax(diff_[j,:]<0)
                
        except:
            out1[j]=float('nan')
            out2[j]=float('nan')
            continue
        try:
            out2[j]=(np.argmax((-diff_[j,int(out1[j]):-1])<(diff_[j,int(out1[j])+1]))+out1[j])*60/1000
        except:
            out2[j]=float('nan')
            continue
        if out1[j]<=21:
            out1[j]=float('nan')
    return out1*60/1000,out2


# %%
def getlinehdif(vpth,vptc,val=2.7):
    out1=np.zeros((9,))
    out2=np.zeros((9,))
    out3=np.zeros((9,))
    for j in range(9):
        diff_=vpth[j,:]-vpth[j,15]
        diff_c=vptc[j,:]-vpth[j,15]
        out1[j]=np.argmax(diff_>val)
        out2[j]=np.argmax(diff_c>val)
        if out2[j]>out1[j]:
            out2[j]=float('nan')
            out3[j]=float('nan')
            continue
        diff_=vpth[j,:]-vptc[j,:]
        diff_[0:20]=0
        try:
            out3[j]=np.argmax(diff_[:]<0)
        except:
            out3[j]=float('nan')
        if out3[j]<=21:
            out3[j]=float('nan')
    return out1*60/1000,out2*60/1000,out3*60/1000


# %%
def getlinepa(vpth,vptc,pah,pac,val=1.5):
    out1=np.zeros((9,))
    out2=np.zeros((9,))
    out3=np.zeros((9,))
    for j in range(9):
        diff_=vpth[j,:]-vpth[j,15]
        out1[j]=np.argmax(diff_>val)
        pahm=pah[j,int(out1[j])]
        out2[j]=np.argmax((pac[j,:]-pahm)<0)-1
        out3[j]=np.argmax((pah[j,:]-pac[j,:])>0)
    return out1*60/1000,out2*60/1000,out3*60/1000


# %%
for i in range(10,20):
    plt.figure()
    diff=pa_s[i,:,:,0]-pa_s[i,:,:,1]
    plt.imshow(diff[:].T,extent=(12,21,0,5),origin='lower',cmap='terrain',vmin=-20,vmax=20)
    plt.colorbar()

# %%
a

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
