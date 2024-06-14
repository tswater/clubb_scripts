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
import os
import rasterio
import seaborn as sns
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
sns.set_theme()
plt.rcParams.update({'figure.max_open_warning': 0})

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
import pickle
fp=open('les_vpt.p','rb')
vpt_s=pickle.load(fp)
vpt_a=pickle.load(fp)
pa_s=pickle.load(fp)
pa_a=pickle.load(fp)
circ=pickle.load(fp)
fp.close()

# %%
vpt_std=np.zeros((len(int_days),9))
j=0
for day in int_days:
    print(day,end=',')
    fdir='/home/tsw35/tyche/data/LES_FULL/fr2_'+str(day)+'_00/'
    hfiles=os.listdir(fdir)
    hfiles.sort()
    i=0
    for file in hfiles:
        if i<5 or i>13:
            i=i+1
            continue
        lall=nc.Dataset(fdir+file,'r')
        vptt=lall['AVV_THV'][0,i,:,:]
        vpt_std[j,i-5]=np.std(vptt)
        #print(str(vpt_std[j,i])[0:4],end=',')
        i=i+1
    j=j+1


# %%
def read_sfc_data(nt,stdt,st=5):
    output=np.zeros((nt,20,20))
    for t in range(nt):
        dt = stdt+datetime.timedelta(hours=t+st)
        tnm = dt.strftime('%Y_%m_%d_%H')
        file = '/home/tsw35/soteria/clubb/data/surfaces_5k/lw/'+tnm+'.tif'
        output[t,:,:] = rasterio.open(file).read(1)
    return output


# %%
ii=0
Tsfc=np.zeros((len(int_days),9,20,20))
for day in int_days:
    print(day,end=',')
    sd=str(day)
    stdt=datetime.datetime(int(sd[0:4]),int(sd[4:6]),int(sd[6:8]),12,0)
    lw=read_sfc_data(9,stdt)
    Tsfc[ii,:]=(lw/(5.67*10**(-8)))**(1/4)
    ii=ii+1

# %%
ii=0
tsfc_std=np.zeros((len(int_days),9))
for day in int_days:
    print(day,end=',')
    tsfc_std[ii,:]=np.std(Tsfc[ii,:].reshape(9,20*20),axis=1)
    ii=ii+1

# %%
plt.hist(tsfc_std.flatten()*2,bins=np.linspace(0,10))

# %%
np.median(tsfc_std.flatten()*2)

# %%
plt.hist(vpt_std.flatten()*2,bins=np.linspace(0,5))

# %%
np.median(vpt_a[:,:,20,1]-vpt_a[:,:,20,0])

# %%
les_dir= '/home/tsw35/tyche/data/LES_1C/'
day=days_nm[0]
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
    plt.figure(figsize=(8,6))
    
    diff=fpht['v2'][1:,0:166]+fpht['u2'][1:,0:166]-fphm['v2'][1:,0:166]-fphm['u2'][1:,0:166]
    avg=np.max(diff,axis=1)
    diff=diff/avg[:,None]
    
    stren=np.mean(fpht['wth'][20:-20,1])
    
    vmax=np.percentile(np.abs(diff),95)
    plt.imshow(diff.T,origin='lower',vmin=0,vmax=1, cmap='terrain',extent=(7,22,0,10))
    plt.xticks(np.arange(7,22,1))
    plt.yticks(np.arange(1,10,1),labels=np.arange(.5,5,.5))
    plt.colorbar(label=r"$\sigma_u^2$ (Circulation)")
    m0,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),2.8)
    m03,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),3)
    m02,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),2)
    m01,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),1)
    m04,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),4)
    m05,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),5)
    m06,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),6)
    m07,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),7)
    m08,maxzl,minz = dT2elev(vpt_s[i,:,:,0],vpt_s[i,:,:,1],np.ones((9,)),8)
    
    m2,maxzl,minz = dT2elev(vpt_a[i,:,:,0],vpt_a[i,:,:,1],tsfc_std[i,:],2)
    m15,maxzl,minz = dT2elev(vpt_a[i,:,:,0],vpt_a[i,:,:,1],tsfc_std[i,:],1.5)
    m175,maxzl,minz = dT2elev(vpt_a[i,:,:,0],vpt_a[i,:,:,1],tsfc_std[i,:],1.75)
    #m5,maxzl,minz = dT2elev(vpt_a[i,:,:,0],vpt_a[i,:,:,1],vpt_std[i,:],5)
    #m7,maxzl,minz = dT2elev(vpt_a[i,:,:,0],vpt_a[i,:,:,1],vpt_std[i,:],7)
    m3,maxzl,minz = dT2elev(vpt_a[i,:,:,0],vpt_a[i,:,:,1],vpt_s[i,:,5,0]-vpt_s[i,:,5,1],2)
    
    plt.plot([12,13,14,15,16,17,18,19,20],m01,'--',color='black',alpha=.75,linewidth=1)
    plt.plot([12,13,14,15,16,17,18,19,20],m02,'--',color='black',alpha=.75,linewidth=1)
    plt.plot([12,13,14,15,16,17,18,19,20],m03,'--',color='black',alpha=.75,linewidth=1)
    plt.plot([12,13,14,15,16,17,18,19,20],m04,'--',color='black',alpha=.75,linewidth=1)
    plt.plot([12,13,14,15,16,17,18,19,20],m05,'--',color='black',alpha=.75,linewidth=1)
    plt.plot([12,13,14,15,16,17,18,19,20],m06,'--',color='black',alpha=.75,linewidth=1)
    plt.plot([12,13,14,15,16,17,18,19,20],m07,'--',color='black',alpha=.75,linewidth=1)
    plt.plot([12,13,14,15,16,17,18,19,20],m08,'--',color='black',alpha=.75,linewidth=1)
    
    plt.plot([12,13,14,15,16,17,18,19,20],m0,'-',linewidth=3)
    plt.plot([12,13,14,15,16,17,18,19,20],m2,'k-',linewidth=3)
    plt.plot([12,13,14,15,16,17,18,19,20],m175,'k--',linewidth=3)
    #plt.plot([12,13,14,15,16,17,18,19,20],m25,'k-',linewidth=3)
    #plt.plot([12,13,14,15,16,17,18,19,20],m15,'k-',linewidth=3)
    plt.plot([12,13,14,15,16,17,18,19,20],m15,'k--',linewidth=3)
    #plt.plot([12,13,14,15,16,17,18,19,20],m5,'r-',linewidth=3)
    #plt.plot([12,13,14,15,16,17,18,19,20],m7,'m-',linewidth=3)
    #plt.plot([12,13,14,15,16,17,18,19,20],maxzl,'k--*')
    plt.plot([12,13,14,15,16,17,18,19,20],minz,'k--o',linewidth=3)
    #plt.legend([r'$z_{hot}$ at $\theta_{v max}$',r'$z_{cold}$ at $\theta_{v max}$',r'$z$ at $\theta_{v crit}$'],framealpha=0.99)
    plt.xlabel('Time (hr)')
    plt.ylabel('Elevation')
    plt.title(str(day)+': '+str(stren)[0:5])
    
    plt.xlim(11.5,20.5)
    
    i=i+1


# %%
def dT2elev(vpth,vptc,vstd,c3):
    out1=np.zeros((9,))
    out2=np.zeros((9,))
    out3=np.zeros((9,))
    for j in range(9):
        diff_=vpth[j,:]-vpth[j,15]
        diff_c=vptc[j,:]-vpth[j,15]
        val=vstd[j]*c3
        #val=(vptc[j,5]-vpth[j,5])*c3
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

# %%
dlht='/home/tsw35/xTyc_shared/clasp/fr2/fr2_20160625_00/'
dlhm='/home/tsw35/xTyc_shared/clasp/fr2/fr2_20160625_01/'
tfile='diag_d01_2016-06-25_20:00:00'

# %%
fpht=nc.Dataset(dlht+tfile,'r')
fphm=nc.Dataset(dlhm+tfile,'r')
alt=fpht['AVP_Z'][0,:]
uht=fpht['AVV_U'][0,:,:,:]
uhm=fphm['AVV_U'][0,:,:,:]
vht=fpht['AVV_V'][0,:,:,:]
vhm=fphm['AVV_V'][0,:,:,:]
wht=fpht['AVV_W'][0,:,:,:]
whm=fphm['AVV_W'][0,:,:,:]

# %%
uTt=np.sqrt(uht**2+vht**2)
uTm=np.sqrt(uhm**2+vhm**2)

# %%
v1=np.mean(np.mean(uTt-uTm,axis=1),axis=1)

# %%
plt.plot(v1[0:166],alt[0:166])

# %%
sigma=10
datav=np.mean(vht[0:5,:,:],axis=0)
data=sci.ndimage.filters.gaussian_filter(datav,sigma,mode='constant')
plt.imshow(data,vmin=-15,vmax=15,origin='lower',cmap='coolwarm')
plt.colorbar()

# %%
data=np.zeros((12,166))
for t in range(12,24):
    tfile='diag_d01_2016-06-25_'+str(t)+':00:00'
    fpht=nc.Dataset(dlhm+tfile,'r')
    uht=fpht['AVV_U'][0,:,:,:]
    data[t-12,:]=np.mean(uht[:,:,:],axis=(1,2))

# %%

# %%
plt.imshow(np.abs(data.T),vmin=0,vmax=10,origin='lower',cmap='terrain',extent=(7,19,0,5))
plt.colorbar()

# %%
plt.plot(data[5,:],alt[0:166])
plt.plot(data[10,:],alt[0:166])

# %%
dlht='/home/tsw35/xTyc_shared/clasp/fr2/fr2_20160625_00/'
dlhm='/home/tsw35/xTyc_shared/clasp/fr2/fr2_20160625_01/'
tfile='diag_d01_2016-06-25_21:00:00'
fpht=nc.Dataset(dlht+tfile,'r')
fphm=nc.Dataset(dlhm+tfile,'r')
alt=fpht['AVP_Z'][0,:]
uht=fpht['AVV_U'][0,:,:,:]
uhm=fphm['AVV_U'][0,:,:,:]
vht=fpht['AVV_V'][0,:,:,:]
vpt=fpht['AVV_THV'][0,:,:,:]

# %%
sigma=10
w_reso=40
zlist=[4,6,8,10,12,16,20,25,30,35,40,45,50,60,70,75,80,85,88,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,115]
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

# %%
### VIDEO ATTEMPT ####
fig = plt.figure(figsize=(6,4.5),dpi=200)
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
cm=plt.cm.get_cmap('Blues')

def animate(i):
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
    ax2.set_ylabel('Elevation ($m$)')
    ax2.set_xlabel('VPT ($K$)')
    ax2.set_xlim(312,320.1)
    ax2.set_ylim(-25,3550)
    fig.subplots_adjust(wspace=.45, hspace=.45)
    
    return fig
    
ani=FuncAnimation(fig,animate,frames=frames,interval=200,repeat=True)
#FFwriter = animation.FFMpegWriter(fps=10)
#ani.save('ani_1.mp4',writer=FFwriter)
HTML(ani.to_jshtml())
#HTML(ani.to_html5_video())

# %%
FFwriter = matplotlib.animation.FFMpegWriter(fps=5)
ani.save('ani_2cols.mp4',writer=FFwriter)

# %%
abmax=np.percentile(np.abs(vpt),90)
for i in range(0,10,1):
    plt.figure(figsize=(4,3))
    sigma=10
    datau=np.mean(vpt[i:i+1,:,:],axis=0)
    data=sci.ndimage.filters.gaussian_filter(datau,sigma,mode='reflect')
    plt.imshow(data,origin='lower',cmap='coolwarm')
    plt.colorbar()
    plt.title('Temperature at '+str((i)*30)+'m')

# %%
plt.show()

# %%
abmax=np.percentile(np.abs(uht),90)
for i in range(0,155,10):
    plt.figure(figsize=(4,3))
    sigma=10
    datav=np.mean(vht[i:i+10,:,:],axis=0)
    datav=sci.ndimage.filters.gaussian_filter(datav-np.mean(datav),sigma,mode='reflect')
    datau=np.mean(uht[i:i+10,:,:],axis=0)
    data=sci.ndimage.filters.gaussian_filter(datau,sigma,mode='reflect')
    plt.imshow(data,origin='lower',cmap='PRGn')
    plt.colorbar()
    plt.title('Velocity at '+str((i+5)*30)+'m')
    break

# %%
plt.imshow(d_vpt[4,:,:],origin='lower',cmap='coolwarm')
x, y = np.meshgrid(np.linspace(0, 500, 10),np.linspace(0, 500, 10))
wind=np.zeros((10,10))
wind_v=np.zeros((10,10))
for i in range(0,500,50):
    for j in range(0,500,50):
        wind[int(i/50),int(j/50)]=data[i,j]
        wind_v[int(i/50),int(j/50)]=datav[i,j]
plt.quiver(x,y,wind,wind_v,color='black',width=.010,alpha=.7)
plt.grid(False)
plt.xticks([])
plt.yticks([])

# %%
plt.quiver(x,y,wind,wind_v)

# %%
plt.imshow(np.mean(wht[30:35,:,:]-whm[30:35,:,:],axis=0),origin='lower',cmap='coolwarm',vmin=-8,vmax=8)
plt.colorbar()

# %%
plt.imshow(fpht['AVV_THV'][0,5,:,:],origin='lower',cmap='coolwarm')
plt.title('Near Surface Temperature FR3')

# %%
plt.imshow(fpht['AVS_TSK'][0,:,:],origin='lower',cmap='coolwarm')
plt.title('Sensible Heat Flux')

# %%
plt.imshow(fpht['AVS_TSK'][0,:,:],origin='lower',cmap='coolwarm')
plt.title('Sensible Heat Flux')

# %%
fr2file='/home/tsw35/tyche/data/LES_FULL/fr2_20160625_00/diag_d01_2016-06-25_20:00:00'
fpfr2=nc.Dataset(fr2file,'r')

# %%
data=fpfr2['AVV_THV'][0,5,0:45,:]
data1=sci.ndimage.filters.gaussian_filter(data,[0,50])
plt.imshow(data1,origin='lower',cmap='coolwarm')
plt.title('Near Surface Temperature FR2')

# %%
plt.figure(figsize=(9,3))
plt.subplot(1,3,1)
plt.imshow(fpfr2['AVS_SH'][0,:,:],origin='lower',cmap='coolwarm')
plt.title('Sensible Heat Flux')
plt.subplot(1,3,2)
plt.imshow(fpht['AVV_THV'][0,5,:,:],origin='lower',cmap='coolwarm')
plt.title('Near Surface\n Temperature FR3')
plt.subplot(1,3,3)
plt.imshow(fpfr2['AVV_THV'][0,5,:,:],origin='lower',cmap='coolwarm')
plt.title('Near Surface\n Temperature FR2')

# %%
j=0
sigma=20
dirlist=os.listdir(lesdir)
dirlist.sort()
dirlist1=[]
print(len(dirlist),flush=True)
for day in dirlist:
    if '_01' in day:
        continue
    if '2017' in day:
        dirlist1.append(day)
lesdir='/home/tsw35/xTyc_shared/clasp/fr2/'
for day in dirlist1:
    try:
        fname=lesdir+day+'/'+'diag_d01_'+day[4:8]+'-'+day[8:10]+'-'+day[10:12]+'_20:00:00'
        fpht=nc.Dataset(fname,'r')
    except Exception as e:
        print(e)
        continue
    uht=fpht['AVV_U'][0,:,:,:]
    vht=fpht['AVV_V'][0,:,:,:]
    H = fpht['AVS_SH'][0,:,:]
    th = fpht['AVV_THV'][0,5,:,:]
    plt.figure(figsize=(12,4))
    plt.subplot(1,3,1)
    data=np.mean(vht[0:5,:,:],axis=0)
    data=sci.ndimage.filters.gaussian_filter(data,sigma)
    plt.imshow(data,origin='lower',cmap='Spectral')
    plt.colorbar()
    plt.subplot(1,3,2)
    data=np.mean(uht[0:5,:,:],axis=0)
    data=sci.ndimage.filters.gaussian_filter(data,sigma)
    plt.imshow(data,origin='lower',cmap='Spectral')
    plt.title(day)
    plt.colorbar()
    plt.subplot(1,3,3)
    plt.imshow(sci.ndimage.filters.gaussian_filter(th,sigma),origin='lower',cmap='seismic')
    plt.title(day)
    plt.colorbar()
    j=j+1

# %%

# %%
j=0
lesdir='/home/tsw35/xTyc_shared/clasp/fr2/'
sigma=25
for day in os.listdir(lesdir):
    try:
        if '2016' in day:
            pass
        else:
            continue
        if '_00' in day:
            fname=lesdir+day+'/'+'diag_d01_'+day[4:8]+'-'+day[8:10]+'-'+day[10:12]+'_20:00:00'
            fpht=nc.Dataset(fname,'r')
        else:
            continue
        day2=day.replace('_00','_01')
        fname2=lesdir+day2+'/'+'diag_d01_'+day[4:8]+'-'+day[8:10]+'-'+day[10:12]+'_20:00:00'
        fphm=nc.Dataset(fname2,'r')
            
    except Exception as e:
        print(e)
        continue
    uht=fpht['AVV_U'][0,:,:,:]
    vht=fpht['AVV_V'][0,:,:,:]
    H = fpht['AVS_SH'][0,:,:]
    th = fpht['AVV_THV'][0,5,:,:]
    dataum=np.mean(fphm['AVV_U'][0,0:5,:,:],axis=0)
    datavm=np.mean(fphm['AVV_V'][0,0:5,:,:],axis=0)
    
    #u2=sci.ndimage.filters.gaussian_filter(,sigma,mode='constant')

    plt.figure(figsize=(12,12))
    plt.subplot(3,3,1)
    datav=np.mean(vht[0:5,:,:],axis=0)
    plt.imshow(datav,vmin=-np.max(np.abs(datav)),vmax=np.max(np.abs(datav)),origin='lower',cmap='Spectral')
    plt.colorbar()
    plt.subplot(3,3,2)
    datau=np.mean(uht[0:5,:,:],axis=0)
    plt.imshow(datau,vmin=-np.max(np.abs(datau)),vmax=np.max(np.abs(datau)),origin='lower',cmap='Spectral')
    plt.title(day)
    plt.colorbar()
    plt.subplot(3,3,3)
    plt.imshow(H,origin='lower',cmap='seismic')
    plt.title(day)
    plt.colorbar()

    plt.subplot(3,3,4)
    data=sci.ndimage.filters.gaussian_filter(datav-np.mean(datavm),sigma,mode='constant')
    plt.imshow(data,vmin=-np.max(np.abs(data)),vmax=np.max(np.abs(data)),origin='lower',cmap='Spectral')
    plt.colorbar()
    plt.subplot(3,3,5)
    data=sci.ndimage.filters.gaussian_filter(datau-np.mean(dataum),sigma,mode='constant')
    plt.imshow(data,vmin=-np.max(np.abs(data)),vmax=np.max(np.abs(data)),origin='lower',cmap='Spectral')
    plt.title(day)
    plt.colorbar()
    
    plt.subplot(3,3,6)
    plt.imshow(th,vmin=np.percentile(th,1),vmax=np.percentile(th,99),origin='lower',cmap='seismic')
    plt.title(day)
    plt.colorbar()
    
    #### HMG ####
    #datav=np.mean(vht[0:5,:,:],axis=0)
    th = fphm['AVV_THV'][0,5,:,:]
    
    plt.subplot(3,3,7)
    data=sci.ndimage.filters.gaussian_filter(datavm-np.mean(datavm),sigma,mode='constant')
    plt.imshow(data,vmin=-np.max(np.abs(data)),vmax=np.max(np.abs(data)),origin='lower',cmap='Spectral')
    plt.colorbar()
    plt.subplot(3,3,8)
    data=sci.ndimage.filters.gaussian_filter(dataum-np.mean(dataum),sigma,mode='constant')
    plt.imshow(data,vmin=-np.max(np.abs(data)),vmax=np.max(np.abs(data)),origin='lower',cmap='Spectral')
    plt.title(day)
    plt.colorbar()
    
    plt.subplot(3,3,9)
    plt.imshow(th,vmin=np.percentile(th,1),vmax=np.percentile(th,99),origin='lower',cmap='seismic')
    plt.title(day)
    plt.colorbar()
    
    if j==10:
        break
    j=j+1

# %%
lesdir='/home/tsw35/xTyc_shared/clasp/fr2/'
break
sigma=20
dirlist=os.listdir(lesdir)
dirlist.sort()
uu={}
dthv={}
vpt={}
dtsk={}
lhet={}
dirlist1=[]
for day in dirlist:
    if '_01' in day:
        continue
    if '20160625' in day:
        dirlist1.append(day)
print(len(dirlist1))
for day in dirlist1:
    print(day,end=':')
    dst=day[4:12]
    dirlist2=os.listdir(lesdir+day)
    dirlist2.sort()
    nt=len(dirlist2)
    uu[dst]=np.zeros((nt,226))
    dthv[dst]=np.zeros((nt,2))
    vpt[dst]=np.zeros((nt,2,226))
    dtsk[dst]=np.zeros((nt,2))
    lhet[dst]=0
    
    #### Calculate Lengthscale of Heterogeneity
    hggt=5
    stdt=datetime.datetime(int(day[4:8]),int(day[8:10]),int(day[10:12]),12)
    Hg,Hv = read_sfc_data('sh',nt,stdt)
    
    Hgg = Hg[hggt,:,:]
    lhet[day] = estimate_l_het(-1,Hgg)
    
    for t in range(nt):
        print(t,end=',')
        try:
            fname=lesdir+day+'/'+dirlist2[t]
            fpht=nc.Dataset(fname,'r')
            day2=day.replace('_00','_01')
            fname2=lesdir+day2+'/'+dirlist2[t]
            fphm=nc.Dataset(fname2,'r')
        except Exception as e:
            print(e)
            continue
        uht=fpht['AVV_U'][0,:,:,:]
        vht=fpht['AVV_V'][0,:,:,:]
        uhm=fphm['AVV_U'][0,:,:,:]
        vhm=fphm['AVV_V'][0,:,:,:]
        #H = fpht['AVS_SH'][0,:,:]
        th = fpht['AVV_THV'][0,:,:,:]
        tsk = fpht['AVS_TSK'][0,:,:]
        
        dtsk[day][t]=2*np.std(tsk)
        
        thvg = np.zeros(th.shape)

        for z in range(uht.shape[0]):
            uhtn=sci.ndimage.filters.gaussian_filter(uht[z,:]-np.mean(uhm[z,:]),sigma)
            vhtn=sci.ndimage.filters.gaussian_filter(vht[z,:]-np.mean(vhm[z,:]),sigma)

            uhmn=sci.ndimage.filters.gaussian_filter(uhm[z,:]-np.mean(uhm[z,:]),sigma,mode='constant')
            vhmn=sci.ndimage.filters.gaussian_filter(vhm[z,:]-np.mean(vhm[z,:]),sigma,mode='constant')

            ubht = np.sqrt(uhtn**2+vhtn**2)
            ubhm = np.sqrt(uhmn**2+vhmn**2)

            uut=np.percentile(ubht,90)
            uum=np.percentile(ubhm,90)

            uu[day][t,z]=uut-uum
            
            thvg[z,:,:]=sci.ndimage.filters.gaussian_filter(th[z,:,:],10)
            
        dthv[day][t]=2*np.std(thvg[10,:,:])
        
        clst=np.zeros(thvg[0,:,:].shape)
        clst[thvg[10,:,:]<np.percentile(thvg[10,:,:],25)]=1
        clst[thvg[10,:,:]>np.percentile(thvg[10,:,:],85)]=2
        if np.all(clst==0) or np.all(clst<=1):
            vpt[day][t,0,:]=np.mean(thvg[0:226,:,:],axis=(1,2))
            vpt[day][t,1,:]=np.mean(thvg[0:226,:,:],axis=(1,2))
        else:
            vpt[day][t,0,:]=np.mean(thvg[0:226,clst==1],axis=1)
            vpt[day][t,1,:]=np.mean(thvg[0:226,clst==2],axis=1)
        
    print()
    alt=fpht['AVP_Z'][0,:]
    odir='/home/tsw35/soteria/clubb/data/les_param/'
    fpo=nc.Dataset(odir+dst+'.nc','w')
    fpo.createDimension('z',size=226)
    fpo.createVariable('z','d',dimensions=('z'))
    fpo['z'][:]=alt[:]
    fpo.createDimension('clst',size=2)
    fpo.createDimension('time',size=nt)
    fpo.createVariable('vpt','d',dimensions=('time','clst','z'))
    fpo.createVariable('uu','d',dimensions=('time','z'))
    fpo.createVariable('dthv','d',dimensions=('time','clst'))
    fpo.createVariable('dtsk','d',dimensions=('time','clst'))
    fpo.createVariable('lhet','d',dimensions=('time'))
    fpo['vpt'][:]=vpt[dst][:]
    fpo['uu'][:]=uu[dst][:]
    fpo['dtsk'][:]=dtsk[dst][:]
    fpo['dthv'][:]=dthv[dst][:]
    fpo['lhet'][:]=lhet[dst][:]
    fpo.close()

# %%

# %%
u_test=uu['fr2_20160625_00'][:]
u2=u_test/np.max(u_test,axis=1)[:,None]

# %%
plt.imshow(u2[:,0:166].T,cmap='terrain',origin='lower',extent=(7,22,0,10))
plt.colorbar()

# %%
a=np.mean(thvg[0:226,clst==1],axis=1)
a.shape

# %%
fpht['AVV_U'][:].shape

# %%
lesdir='/home/tsw35/xTyc_shared/clasp/fr2/'
ZZ=np.meshgrid(np.linspace(0,519,520),np.linspace(0,519,520),np.linspace(0,165,166))
x=ZZ[0].flatten()
msk00=np.random.randint(0,len(x),1000000)
dataht=np.zeros((16,1000000))
datahm=np.zeros((16,1000000))
Hht=np.zeros((16,520,520))
Hhm=np.zeros((16,520,520))
x=x[msk00]
y=ZZ[1].flatten()[msk00]
z=ZZ[2].flatten()[msk00]

anidir1=lesdir+'fr2_20170629_01'
anidir0=lesdir+'fr2_20170629_00'
ani1list=os.listdir(anidir1)
ani1list.sort()
ani0list=os.listdir(anidir0)
ani0list.sort()
i=0
for file in ani0list:
    print('.',end='')
    fp=nc.Dataset(anidir0+'/'+file,'r')
    dataht[i,:]=np.transpose(fp['AVV_QC'][0,0:166,:,:],[1,2,0]).flatten()[msk00]
    Hht[i,:]=fp['AVS_SH'][0,:,:]
    i=i+1
i=0
for file in ani1list:
    print('*',end='')
    fp=nc.Dataset(anidir1+'/'+file,'r')
    datahm[i,:]=np.transpose(fp['AVV_QC'][0,0:166,:,:],[1,2,0]).flatten()[msk00]
    Hhm[i,:]=fp['AVS_SH'][0,:,:]
    i=i+1

# %%
vmax=np.percentile(dataht[dataht>0],99)
zsh=np.ones((520,520))
ZZ2=np.meshgrid(np.linspace(0,519,520),np.linspace(0,519,520))

minn, maxx = Hht.min(), np.percentile(Hht,99)
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='coolwarm')
m.set_array([])

# %%
### VIDEO ATTEMPT ####
fig = plt.figure(figsize=(8,12),facecolor='silver')
ax1=fig.add_subplot(211,projection='3d')
ax2=fig.add_subplot(212,projection='3d')
cm=plt.cm.get_cmap('Blues')

def animate(i):
    ax1.clear()
    qc=dataht[i,:]
    xx=x[qc>0]
    yy=y[qc>0]
    zz=z[qc>0]
    qc=qc[qc>0]
    
    H=Hht[i,:]
    fcolors = m.to_rgba(H)
    ax1.plot_surface(ZZ2[0],ZZ2[1],zsh,rstride=10,cstride=10,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)
    
    sc=ax1.scatter(xx,yy,zz,s=10,alpha=.15,c=qc,cmap=cm,vmin=0,vmax=vmax)
    ax1.set_xlim(0,520*.9)
    ax1.set_ylim(0,520*.9)
    ax1.set_zlim(15,140*.5)
    ax1.set_facecolor("silver")
    
    ax1.view_init(elev=10,azim=-60)
    ax1.axis('off')
    ax1.set_title('Heterogeneous',fontsize=20)
    
    #plt.colorbar(sc)
    
    #if i==0:
    #    fig.colorbar(sc,ax=ax1)
    
    ax2.clear()
    qc=datahm[i,:]
    xx=x[qc>0]
    yy=y[qc>0]
    zz=z[qc>0]
    qc=qc[qc>0]
    
    H=Hhm[i,:]
    fcolors = m.to_rgba(H)
    ax2.plot_surface(ZZ2[0],ZZ2[1],zsh,rstride=10,cstride=10,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)
    
    sc=ax2.scatter(xx,yy,zz,s=5,alpha=.2,c=qc,cmap=cm,vmin=0,vmax=vmax)
    ax2.set_xlim(0,520*.9)
    ax2.set_ylim(0,520*.9)
    ax2.set_zlim(15,140*.5)
    ax2.set_facecolor("silver")
    ax2.view_init(elev=10,azim=-60)
    ax2.axis('off')
    ax2.set_title('Homogeneous',fontsize=20)
    
    #plt.colorbar(sc)
    
    #if i==0:
    #    fig.colorbar(sc,ax=ax2)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    #plt.tight_layout()
    
    return fig
    
ani=FuncAnimation(fig,animate,frames=16,interval=250,repeat=True)
#FFwriter = animation.FFMpegWriter(fps=10)
#ani.save('ani_1.mp4',writer=FFwriter)
HTML(ani.to_jshtml())
#HTML(ani.to_html5_video())

# %%
FFwriter = matplotlib.animation.FFMpegWriter(fps=4)
ani.save('ani_les20170629.mp4',writer=FFwriter)

# %%
fig = plt.figure(figsize=plt.figaspect(1),facecolor='silver',layout="constrained")
subfigs = fig.subfigures(2, 1, hspace=0,frameon=False,)
ax1=subfigs[0].add_subplot(111,projection='3d')
ax2=subfigs[1].add_subplot(111,projection='3d')
cm=plt.cm.get_cmap('Blues')

def animate(i):
    ax1.clear()
    qc=dataht[i,:]
    xx=x[qc>0]
    yy=y[qc>0]
    zz=z[qc>0]
    qc=qc[qc>0]
    
    H=Hht[i,:]
    fcolors = m.to_rgba(H)
    ax1.plot_surface(ZZ2[0],ZZ2[1],zsh,rstride=10,cstride=10,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)
    
    sc=ax1.scatter(xx,yy,zz,s=10,alpha=.15,c=qc,cmap=cm,vmin=0,vmax=vmax)
    ax1.set_xlim(0,520)
    ax1.set_ylim(0,520)
    ax1.set_zlim(0,140)
    ax1.set_facecolor("blue")
    
    ax1.view_init(elev=10,azim=-60)
    ax1.axis('off')
    ax1.set_title('Heterogeneous',fontsize=20)
    #ax1.set_box_aspect((52,52,14))
    
    #plt.colorbar(sc)
    
    #if i==0:
    #    fig.colorbar(sc,ax=ax1)
    
    ax2.clear()
    qc=datahm[i,:]
    xx=x[qc>0]
    yy=y[qc>0]
    zz=z[qc>0]
    qc=qc[qc>0]
    
    H=Hhm[i,:]
    fcolors = m.to_rgba(H)
    ax2.plot_surface(ZZ2[0],ZZ2[1],zsh,rstride=10,cstride=10,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)
    
    sc=ax2.scatter(xx,yy,zz,s=5,alpha=.2,c=qc,cmap=cm,vmin=0,vmax=vmax)
    ax2.set_xlim(0,520)
    ax2.set_ylim(0,520)
    ax2.set_zlim(30,60)
    ax2.set_facecolor("blue")
    ax2.view_init(elev=10,azim=-60)
    ax2.axis('off')
    ax2.set_title('Homogeneous',fontsize=20)
    rt=1
    ax2.set_box_aspect((52*rt,52*rt,14))
    #ax2.set_box_aspect((1,1,1))
    
    
    #plt.colorbar(sc)
    
    #if i==0:
    #    fig.colorbar(sc,ax=ax2)
    #fig.subplots_adjust(left=-.5, right=2, bottom=0, top=1)
    #fig.subplots_adjust(wspace=.1, hspace=0)
    #plt.tight_layout()
    
    
    return fig

animate(12)

# %%

# %%

# %%
from mpl_toolkits.mplot3d import Axes3D
import matplotlib

# %%
qc=fpht['AVV_QC'][0,:,:,:]
qcT=np.transpose(qc,[1,2,0])
#fig = plt.figure()
#ax=fig.add_subplot(projection='3d')
ZZ=np.meshgrid(np.linspace(0,519,520),np.linspace(0,519,520),np.linspace(0,225,226))

# %%
qc=fpht['AVV_QC'][0,0:166,:,:]
qcT=np.transpose(qc,[1,2,0])
ZZ=np.meshgrid(np.linspace(0,519,520),np.linspace(0,519,520),np.linspace(0,165,166))
fig = plt.figure(figsize=(12,6))
ax=fig.add_subplot(111,projection='3d')
cm=plt.cm.get_cmap('Blues')
qcc=qcT.flatten()
qc2=qcc[qcc>0]
msk=np.random.randint(0,len(qc2),100000)
H=fpht['AVS_SH'][0,:,:]
zsh=np.ones((520,520))
ZZ2=np.meshgrid(np.linspace(0,519,520),np.linspace(0,519,520))

minn, maxx = H.min(), H.max()
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='coolwarm')
m.set_array([])
fcolors = m.to_rgba(H)

#ax.plot_surface(ZZ2[0],ZZ2[1],zsh,rstride=10,cstride=10,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)
ax.plot_surface(ZZ2[0],ZZ2[1],zsh,rstride=10,cstride=10,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)

vmax=np.percentile(qc2,90)

x=ZZ[0].flatten()[qcc>0][msk]
y=ZZ[1].flatten()[qcc>0][msk]
z=ZZ[2].flatten()[qcc>0][msk]

#ax.grid(False)

sc=ax.scatter(x,y,z,s=10,alpha=.25,c=qc2[msk],cmap=cm,vmax=vmax)
ax.view_init(elev=15,azim=-60)
ax.axis('off')

# %%
ax.azim

# %%
# sci.ndimage.filters.gaussian_filter(datavm-np.mean(datavm),sigma,mode='constant')
lesdir='/home/tsw35/xTyc_shared/clasp/fr2/'
sigma=10
vpt=[]
j=0
for day in os.listdir(lesdir):
    try:
        if '2016' in day:
            pass
        else:
            continue
        if '_00' in day:
            fname=lesdir+day+'/'+'diag_d01_'+day[4:8]+'-'+day[8:10]+'-'+day[10:12]+'_20:00:00'
            fpht=nc.Dataset(fname,'r')
        else:
            continue
        day2=day.replace('_00','_01')
        #fname2=lesdir+day2+'/'+'diag_d01_'+day[4:8]+'-'+day[8:10]+'-'+day[10:12]+'_20:00:00'
        #fphm=nc.Dataset(fname2,'r')
            
    except Exception as e:
        print(e)
        continue
    uht=fpht['AVV_U'][0,:,:,:]
    vht=fpht['AVV_V'][0,:,:,:]
    #uhm=fphm['AVV_U'][0,:,:,:]
    #vhm=fphm['AVV_V'][0,:,:,:]
    
    #H = fpht['AVS_SH'][0,:,:]
    plt.figure()
    plt.plot(np.mean(uht[0:166,:,:],axis=(1,2)),alt[0:166])
    plt.plot(np.mean(vht[0:166,:,:],axis=(1,2)),alt[0:166])
    plt.title(day)

# %%
plt.plot(np.mean(uht[0:166,:,:],axis=(1,2)),alt)
plt.ylim(0,5000)
plt.xlim(-2,2)

# %%
# sci.ndimage.filters.gaussian_filter(datavm-np.mean(datavm),sigma,mode='constant')
lesdir='/home/tsw35/xTyc_shared/clasp/fr2/'
sigma=10
vpt=[]
j=0
for day in os.listdir(lesdir):
    try:
        if '20160625' in day:
            pass
        else:
            continue
        if '_00' in day:
            fname=lesdir+day+'/'+'diag_d01_'+day[4:8]+'-'+day[8:10]+'-'+day[10:12]+'_20:00:00'
            fpht=nc.Dataset(fname,'r')
        else:
            continue
        day2=day.replace('_00','_01')
        fname2=lesdir+day2+'/'+'diag_d01_'+day[4:8]+'-'+day[8:10]+'-'+day[10:12]+'_20:00:00'
        fphm=nc.Dataset(fname2,'r')
            
    except Exception as e:
        print(e)
        continue
    #uht=fpht['AVV_U'][0,:,:,:]
    #vht=fpht['AVV_V'][0,:,:,:]
    #H = fpht['AVS_SH'][0,:,:]
    thv = fpht['AVV_THV'][0,0:166,:,:]
    thvg = np.zeros(thv.shape)
    vpt.append(np.zeros((2,166)))
    
    for z in range(166):
        thvg[z,:,:]=sci.ndimage.filters.gaussian_filter(thv[z,:,:],sigma)
        #print('.',end='')
    clst=np.zeros(thvg[0,:,:].shape)
    clst[thvg[10,:,:]<np.percentile(thvg[10,:,:],20)]=1
    clst[thvg[10,:,:]>np.percentile(thvg[10,:,:],80)]=2
    for z in range(166):
        vpt[j][0,z]=np.mean(thvg[z,clst==1])
        vpt[j][1,z]=np.mean(thvg[z,clst==2])
    print(j)
    j=j+1
vpt=np.array(vpt)

# %%
u1=[-2,-4]
u2=[2,4]
np.dot(u2,u1)/np.sqrt(u1[0]**2+u1[1]**2)

# %%
np.sqrt(u2[0]**2+u2[1]**2)

# %%
a=.8
b=.9
c=.2
d=.1
aa=np.sqrt(a**2+b**2)-np.sqrt(c**2+d**2)
bb=np.sqrt((a-c)**2+(b-d)**2)
print(aa)
print(bb)

# %%
vpt=vpto.copy()

# %%
j=6
for i in range(vpt.shape[0]):
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    plt.subplot(1,6,j+1)
    plt.plot(vpt[i,0,:],alt[0:166])
    plt.plot(vpt[i,1,:],alt[0:166])
    plt.title(days2017[i])
    plt.xlim(np.min(vpt[i,0,:])-1,np.min(vpt[i,0,:])+7)
    j=j+1

# %%
plt.imshow(thvg[110,:,:])
plt.colorbar()
plt.figure()
plt.imshow(clst)

# %%
days2017=[]
for day in os.listdir(lesdir):
    if '2017' in day:
        pass
    else:
        continue
    if '_00' in day:
        days2017.append(day[4:16])

# %%
vpt=np.zeros((2,166))
clst=np.zeros(thvg[0,:,:].shape)
clst[thvg[10,:,:]<np.percentile(thvg[10,:,:],25)]=1
clst[thvg[10,:,:]>np.percentile(thvg[10,:,:],75)]=2
for z in range(166):
    vpt[0,z]=np.mean(thvg[z,clst==1])
    vpt[1,z]=np.mean(thvg[z,clst==2])

# %%
plt.figure(figsize=(3,6))
plt.plot(vpt[0,:],alt[0:166])
plt.plot(vpt[1,:],alt[0:166])

# %%

# %%
sfc_dir   = '/home/tsw35/soteria/clubb/data/surfaces_5k/'
dx=5000

#### CALCULATE LENGTHSCALE OF HETEROGENEITY ####
def estimate_l_het(l_het,Hg_,cut=.05,samp=100000):
    '''
        Parameters

        l_het  : -1: compute full lengthscale, -2: compute parallel to wind
        cut    : percentage cuttoff for lengthscale
        Hg     : a N by N grid of the sensible heat fluxes
        samp   : how many points to sample when constructing cov matrix
                 only applies to -1
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

        idx = np.random.choice(len(H_flat),size=400,replace=False)
        H_sb = H_flat[idx]
        r_Hsb = r_flat[idx]
        c_Hsb = c_flat[idx]
        mu=np.mean(H_sb)
        a = (H_sb[:,np.newaxis].T - mu)*(H_sb[:,np.newaxis]-mu)
        h = ((r_Hsb[:,np.newaxis] - r_Hsb.T)**2 + \
            (c_Hsb[:,np.newaxis] - c_Hsb.T)**2)**0.5
        Qf = a.flatten()
        hf = h.flatten()

        bins = np.linspace(2500,75000-2500,15)
        means=np.zeros((len(bins)-1,))
        for i in range(len(bins)-1):
            means[i]=np.mean(Qf[(hf>bins[i])&(hf<bins[i+1])])

        l_het_ = np.mean(bins[0:-1][means<=(cut*means[0])][0:1])

        #print('Non Dimensional Lengthscale of Heterogeneity')
        #print(str(l_het_))
        #print()
    return l_het_
        
#### READ IN SURFACE DATA ####
def read_sfc_data(var,nt,stdt,override='X'):
    nx=20
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
dir2c='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(dir2c)
filelist.sort()
for file in filelist:
    fp=nc.Dataset(dir2c+file,'r')
    u90=fp['u90'][:]
    u10=fp['u10'][:]
    v90=fp['v90'][:]
    v10=fp['v10'][:]
    
    plt.figure(figsize=(15,3))
    plt.subplot(1,4,1)
    plt.title('u90: '+str(file))
    plt.imshow(u90[:,0:166].T,cmap='terrain',origin='lower',extent=(7,22,0,10))
    plt.colorbar()
    plt.subplot(1,4,2)
    plt.title('u10: '+str(file))
    plt.imshow(u10[:,0:166].T,cmap='terrain',origin='lower',extent=(7,22,0,10))
    plt.colorbar()
    plt.subplot(1,4,3)
    plt.title('v90: '+str(file))
    plt.imshow(v90[:,0:166].T,cmap='terrain',origin='lower',extent=(7,22,0,10))
    plt.colorbar()
    plt.subplot(1,4,4)
    plt.title('v10: '+str(file))
    plt.imshow(v10[:,0:166].T,cmap='terrain',origin='lower',extent=(7,22,0,10))
    plt.colorbar()
    break

# %%
dir2c='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(dir2c)
filelist.sort()
for file in filelist:
    fp=nc.Dataset(dir2c+file,'r')
    u90=fp['u90'][:]
    u10=fp['u10'][:]
    v90=fp['v90'][:]
    v10=fp['v10'][:]
    
    plt.figure(figsize=(8,3))
    plt.subplot(1,2,1)
    plt.title('u: '+str(file))
    data=u90[:,0:166].T-u10[:,0:166].T
    nch=(u90[:,0:166].T/np.abs(u90[:,0:166].T))==(u10[:,0:166].T/np.abs(u10[:,0:166].T))
    data[nch]=float('nan')
    plt.imshow(data,cmap='coolwarm',origin='lower',extent=(7,22,0,10))
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.title('v: '+str(file))
    data=v90[:,0:166].T-v10[:,0:166].T
    nch=(v90[:,0:166].T/np.abs(v90[:,0:166].T))==(v10[:,0:166].T/np.abs(v10[:,0:166].T))
    data[nch]=float('nan')
    plt.imshow(data,cmap='coolwarm',origin='lower',extent=(7,22,0,10))
    plt.colorbar()

# %%
dir2c='/home/tsw35/soteria/clubb/data/les_param2/'
filelist=os.listdir(dir2c)
filelist.sort()
u1s=[]
u2s=[]
u3s=[]
ums=[]
um2s=[]
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
    fp=nc.Dataset(dir2c+file,'r')
    t1=3
    t2=12
    ht=33
    ht0=5
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
    uu=(u90-u10)/2
    m2=(np.abs(u10)>np.abs(u90))
    u90[m2]=u10[m2]
    u90=np.abs(u90)
    
    msk=((v90/np.abs(v90))==(v10/np.abs(v10)))|(dvpt<0)
    v90[msk]=0
    v10[msk]=0
    meanv=np.abs((v90+v10)/2)
    vv=(v90-v10)/2
    m2=(np.abs(v10)>np.abs(v90))
    v90[m2]=v10[m2]
    v90=np.abs(v90)
    
    u1=np.sqrt(uu**2+vv**2)
    u2=np.sqrt(v90**2+u90**2)
    u3=np.sqrt((v90+meanv)**2+(u90+meanu)**2)
    dvpt=dvpt[u2>0]
    u3=u3[u2>0]
    u2=u2[u2>0]
    u1=u1[u1>0]
    #u2=u2[u2>0]*300/(lhet**(1/2)*9.8**(1/2))
    u2s.extend(u2)
    u1s.extend(u1)
    u3s.extend(u3)
    um=dvpt/300*(lhet**(1/2)*9.8**(1/2))
    um2=dvpt**1.25/300*(lhet**(1/2)*9.8**(1/2))
    ums.extend(um)
    um2s.extend(um2)
    dvpts.extend(dvpt)
    
    u2=np.array(u2)
    u1=np.array(u1)
    um=np.array(um)
    u3=np.array(u3)
    #u22=np.mean(u2,0)
    
    if len(u2)<5:
        continue
    plt.figure()
    plt.scatter(um,u2)
    plt.xlim(0,np.max([um,u2]))
    plt.ylim(0,np.max([um,u2]))
    
    model=LinearRegression(fit_intercept=False)
    model.fit(um.reshape(-1, 1),u2.reshape(-1, 1))
    cf=model.coef_
    sc=model.score(um.reshape(-1, 1),u2.reshape(-1, 1))
    
    plt.title(file+': '+str(cf)+' : '+str(sc))
        
um2s=np.array(um2s)
u3s=np.array(u3s)
dvpts=np.array(dvpts)
u1s=np.array(u1s)
ums=np.array(ums)
u2s=np.array(u2s)

# %%
stats.pearsonr(ums,u2s)[0]

# %%
cs=np.array(u2s)/np.array(ums)
plt.hist(cs,bins=np.linspace(0,5))
print(np.median(cs))

# %%
stats.pearsonr(ums,u2s)[0]

# %%
u1s=np.array(u1s)
ums=np.array(ums)
u2s=np.array(u2s)

# %%
pr1=[]
pr2=[]
pr3=[]
ass=np.linspace(.001,3,100)
for j in ass:
    msss=dvpts>=0
    um3=dvpts[msss]**(j)
    pr3.append(stats.pearsonr(um3,u3s[msss])[0])
    pr2.append(stats.pearsonr(um3,u2s[msss])[0])
    pr1.append(stats.pearsonr(um3,u1s[msss])[0])
plt.plot(ass,pr1)
plt.plot(ass,pr2)
plt.plot(ass,pr3)

# %%
ums=np.array(ums)

# %%
um3=ums/dvpts*dvpts**(.25)

# %%
mssk2=(dvpts>=0)

# %%
u1s=np.array(u1s)

# %%
stats.pearsonr(um3[mssk2],u1s[mssk2])[0]

# %%
import sklearn.metrics
import scipy.stats as stats

# %%
stats.pearsonr(ums,u2s)

# %%
stats.spearmanr(dvpts,u1s)

# %%

# %%
plt.figure(figsize=(3,3))
plt.hexbin(ums*1.2,u2s,gridsize=40,extent=(0,3,0,3),cmap='terrain')
plt.plot([0,3],[0,3],'k--',alpha=.5)
#plt.xlim(0,3)
#plt.ylim(0,3)
plt.xlabel('Model Predicted Velocity')
plt.ylabel('LES Velocity')

# %%
plt.hexbin(um3*a+b,u2s,gridsize=50,cmap='terrain',mincnt=1)
plt.plot([0,5],[0,5],'k--',alpha=.5)

# %%
a,_,_,_ = np.linalg.lstsq(ums[:,np.newaxis], u2s)

# %%
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
model=LinearRegression(fit_intercept=False)
msk=(ums<=20)&(u2s<=20)
model.fit(ums[msk].reshape(-1, 1),u2s[msk].reshape(-1, 1))
print(model.coef_)
print(model.score(ums[msk].reshape(-1, 1),u2s[msk].reshape(-1, 1)))
print(model.intercept_)

# %%
plt.figure(figsize=(4,3),dpi=300)
plt.plot(np.linspace(0,7),np.linspace(0,7),'k--',zorder=0)
plt.scatter(ums[msk]*1.36,u2s[msk],s=6,alpha=.1,zorder=1)
plt.xlim(0,7)
plt.ylim(0,7)
plt.xlabel('Model Velocity ($ms^{-1}$)',fontsize=14)
plt.ylabel('LES Velocity ($ms^{-1}$)',fontsize=14)

# %%
plt.scatter(ums,u2s,s=5,alpha=.25)
#plt.xlim(0,5)
#plt.ylim(0,5)
plt.plot(np.linspace(0,6),np.linspace(0,6)*1.1,'k--')
plt.plot(np.linspace(0,6),np.linspace(0,6)*1.1,'k--')

# %%
ums=np.array(ums)

# %%
plt.hexbin(ums*.75,u1s,gridsize=50,extent=(0,3,0,3),cmap='terrain',mincnt=1)
plt.plot([0,3],[0,3],'k--',alpha=.5)

# %%
plt.figure(figsize=(5,5))
ums=np.array(ums)
#plt.plot([0,5],[0,5],'k--',alpha=.5)
plt.scatter(ums,u2s,s=10,alpha=.25)
plt.plot(np.linspace(0,5),np.linspace(0,5),'k--')
plt.xlim(0,5)
plt.ylim(0,5)

# %%

# %%
plt.hexbin(ums*.6,u2s,gridsize=50,extent=(0,5,0,5),cmap='terrain',mincnt=1)
plt.plot([0,5],[0,5],'k--',alpha=.5)

# %%
cs=np.array(u2s)/np.array(ums)
plt.hist(cs,bins=np.linspace(0,5))
print(np.mean(cs))

# %%
plt.hist(dvpts,bins=np.linspace(0,4))
print()


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
tskdir='/home/tsw35/soteria/clubb/data/surfaces_5k/lw'
tsk={}
for file in filelist:
    stdt=datetime.datetime(int(file[0:4]),int(file[4:6]),int(file[6:8]),12,0)
    lwg,lwv=read_sfc_data('lw',16,stdt)
    tskg=(lwg/(5.67*10**(-8)))**(1/4)
    tsk[file]=np.array([np.std(tskg,axis=(1,2)),np.std(tskg,axis=(1,2))])
    print('.',end='',flush=True)

# %%
fp['dtsk'][:].shape
['sgp_20150801', 'sgp_20160610', 'sgp_20160625', 'sgp_20160716', 'sgp_20170716', 'sgp_20170717', 'sgp_20170802', 'sgp_20180707', 'sgp_20180709', 'sgp_20190707']



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
dir2c='/home/tsw35/soteria/clubb/data/les_param/'
filelist=os.listdir(dir2c)
filelist.sort()
cc1=.85
dayst={}
filelist=['20160625.nc','20170716.nc','20170717.nc']
filelist=['20160625.nc','20170716.nc','20170717.nc','20190707.nc','20180709.nc','20180707.nc','20170802.nc','20150801.nc','20160610.nc','20160716.nc']
for day in filelist:
    dayst[day]=day[0:4]+'-'+day[4:6]+'-'+day[6:8]
for file in filelist:
    fp=nc.Dataset(dir2c+file,'r')
    u90=fp['u90'][:]
    u10=fp['u10'][:]
    v90=fp['v90'][:]
    v10=fp['v10'][:]
    
    vpt=fp['vpt'][:]
    dtsk=2*tsk[file].T
    
    
    maxzh,maxzl,minzh,min0 = circ_h(vpt,dtsk,cc1)
    
    plt.figure(figsize=(4.5,3),dpi=300)
    #plt.subplot(1,2,1)
    plt.title(dayst[file])
    #data=u90[:,0:166].T-u10[:,0:166].T
    data=np.sqrt((u90[:,0:166].T-u10[:,0:166].T)**2+(v90[:,0:166].T-v10[:,0:166].T)**2)
    
    data=data/np.max(data,axis=0)
    
    plt.imshow(data,cmap='terrain',origin='lower',extent=(7,22,0,10))
    plt.colorbar(label=r'Normalized $u_r$')
    
    plt.plot(np.linspace(7,22,16),maxzh*30*2/1000,'k--')
    plt.plot(np.linspace(7,22,16),maxzl*30*2/1000,'k-')
    plt.plot(np.linspace(7,22,16),minzh*30*2/1000,'k:')
    plt.plot(np.linspace(7,22,16),min0*30*2/1000,':',c='grey')
    plt.ylim(0,10)
    plt.xlim(8,20)
    plt.yticks([0,2,4,6,8,10],labels=[0,1,2,3,4,5])
    plt.xticks([8,12,16,20],labels=['8:00','12:00','16:00','20:00'])
    plt.ylabel('Elevation ($km$)')
    plt.xlabel('Local Time')
    plt.legend([r'$z_{max_1}$',r'$z_{max_2}$',r'$z_{crit}$'],loc='upper left',framealpha=.9)
    '''
    plt.subplot(1,2,2)
    plt.title('v: '+str(file))
    data=v90[:,0:166].T-v10[:,0:166].T
    #nch=(v90[:,0:166].T/np.abs(v90[:,0:166].T))==(v10[:,0:166].T/np.abs(v10[:,0:166].T))
    #data[nch]=float('nan')
    plt.imshow(data,cmap='coolwarm',origin='lower',extent=(7,22,0,10))
    plt.colorbar()
    
    plt.plot(np.linspace(7,22,16),maxzh*40*2/1000,'k--')
    plt.plot(np.linspace(7,22,16),maxzl*40*2/1000,'k-')
    plt.plot(np.linspace(7,22,16),minzh*40*2/1000,'k:')
    '''
plt.show()


# %%

# %%
def circ_h(vpt_,dtsk_,cc1=.85):
    maxzh=np.zeros((16,))
    minzh=np.zeros((16,))
    min0=np.zeros((16,))
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
            min0[t]=np.where((vpt_[t,1,:]-vpt_[t,0,:])<0)[0][0]-1
            minzh[t]=np.where((vpt_[t,1,:]-vpt_[t,0,:])<0.25)[0][0]-1
        except:
            minzh[t]=0
            min0[t]=0
    return maxzh,maxzl,minzh,min0


# %%
for file in filelist:
    fp=nc.Dataset(dir2c+file,'r')
    
    vpt=fp['vpt'][:]
    dtsk=fp['dtsk'][:]
    
    data=vpt[:,1,0:166]-vpt[:,0,0:166]
    #data[data>=0]=float('nan')
    plt.figure()
    plt.imshow(data.T,vmin=-np.max(np.abs(data)),vmax=np.max(np.abs(data)),cmap='coolwarm',origin='lower',extent=(7,22,0,10))
    plt.colorbar()
    plt.title(file)
    break

# %%
test=data[10,:]
a=np.where(test<10)[0][0]

# %%
a

# %%

# %%
data.shape

# %%

# %%
