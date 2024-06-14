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
import rasterio
import os
import datetime

# %%
clubb_dir= '/home/tsw35/tyche/clubb/'

# %%
nx=20
dx=5000
stdate  = '2017-07-16T10:00:00.000'# start date iso format
enddate = '2017-07-17T03:00:00.000'# end date iso format
stdt = datetime.datetime.fromisoformat(stdate)
endt = datetime.datetime.fromisoformat(enddate)
nt     = int((endt-stdt).seconds/3600)+1
sfc_dir   = '/home/tsw35/soteria/clubb/data/surfaces_5k/'
k=2


# %%
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


# %%
def compute_width(clst_2d):
    W = np.zeros((k,k))
    cdirect = np.zeros((k,k,2))
    for i in range(nx):
        for j in range(nx):
            cl = int(clst_2d[i,j])
            if i>0:
                cl_2 = int(clst_2d[i-1,j])
                W[cl,cl_2]=W[cl,cl_2]+dx
                cdirect[cl,cl_2,1]=cdirect[cl,cl_2,1]+1
            if i<(nx-1):
                cl_2 = int(clst_2d[i+1,j])
                W[cl,cl_2]=W[cl,cl_2]+dx
                cdirect[cl,cl_2,1]=cdirect[cl,cl_2,1]-1
            if j>0:
                cl_2 = int(clst_2d[i,j-1])
                W[cl,cl_2]=W[cl,cl_2]+dx
                cdirect[cl,cl_2,0]=cdirect[cl,cl_2,0]-1
            if j<(nx-1):
                cl_2 = int(clst_2d[i,j+1])
                W[cl,cl_2]=W[cl,cl_2]+dx
                cdirect[cl,cl_2,0]=cdirect[cl,cl_2,0]+1
    return W,cdirect


# %%
#### COMPUTE DIRECTION OF EACH CONNECTION ####
def circdir(cgrid):
    Wd = np.zeros((2,2))
    tmp=np.zeros((4,))
    # col connections [east west flow]
    for i in range(nx):
        for j in range(nx-1):
            c_1=cgrid[i,j]
            c_2=cgrid[i,j+1]
            if c_1==c_2:
                continue
            elif c_1==0:
                tmp[0]=tmp[0]+1
            else:
                tmp[2]=tmp[2]+1
    # row connections [north south flow]
    for j in range(nx):
        for i in range(nx-1):
            c_1=cgrid[i,j]
            c_2=cgrid[i+1,j]
            if c_1==c_2:
                continue
            elif c_1==0:
                tmp[3]=tmp[3]+1
            else:
                tmp[1]=tmp[1]+1
    Wd[:]=np.sum(tmp)
    circ_d=np.array([tmp[0]+tmp[2],tmp[1]+tmp[3]])
    return Wd,circ_d


# %%
def thresh(Hv_,typ):
    if typ=='median':
        th=np.median(Hv_)
    elif typ=='mean':
        th=np.mean(Hv_)
    elif typ=='test':
        ctff = 0
        prev=0
        diff=100
        idd=0
        for i in range(25,76,2):
            ctf = np.percentile(Hv_,i)
            if (ctf-prev)<diff:
                diff=ctf-prev
                ctff=ctf
                idd=i
            prev=ctf
        th=ctff
        print(idd,flush=True)
    else:
        th=np.percentile(Hv_,typ)
    clst=np.zeros(Hv_.shape)
    clst[Hv_>=th]=0
    clst[Hv_<th]=1
    clst_2=np.reshape(clst,(1,nx,nx))
    return clst,clst_2


# %%
def area(clst):
    A=np.zeros((k))
    A[1]=np.sum(clst)
    A[0]=np.sum(clst==0)
    return A


# %%
Hg,Hv = read_sfc_data('sh',nt,stdt)

# %%
plt.rcParams.update({'figure.max_open_warning': 0})

# %%
clst,clst_2d=thresh(Hv[10,:],'mean')
W,cdirect=compute_width(clst_2d)
A=area(clst)

# %%
print(A)
print(A/W[0,1])

# %%
for i in np.linspace(25,75,21):
    print(i)

# %%
days = ['20180805','20150829','20160719','20170716','20170705','20180916','20190804','20170717','20160619','20180710','20170923']
days.sort()

# %%
day_dir  = '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
days=[]
for file in os.listdir(day_dir):
    if '_01' in file:
        days.append(file[8:16])
    else:
        continue
days.sort()

# %%
i=0
t=7
n=len(days)
A_s  =np.zeros((n,k))
W_s  =np.zeros((n,))
H_s  =np.zeros((n,nx,nx))
clsts=np.zeros((n,nx,nx))
cd_s =np.zeros((n,2))
for day in days:
    file = 'trimdoz_'+day+'_01.nc'
    t_init  = 39600
    t_final = 97200
    pre_dt = datetime.datetime(int(file[8:12]),int(file[12:14]),int(file[14:16]),0,0)
    stdt = pre_dt+datetime.timedelta(seconds=t_init)
    endt = pre_dt+datetime.timedelta(seconds=t_final)
    Hg,Hv = read_sfc_data('sh',nt,stdt)
    H_s[i,:]=Hg[t,:]
    clst,clsts[i,:]=thresh(Hv[t,:],'test')
    W,cdirect=circdir(clsts[i,:])
    W_s[i]=W[0,1]
    cd_s[i,:]=cdirect[:]
    A_s[i,:]=area(clst)*dx*dx
    i=i+1

# %%
i=0
t=7
n=len(days)
A_s  =np.zeros((n,k))
W_s  =np.zeros((n,))
H_s  =np.zeros((n,nx,nx))
clsts=np.zeros((n,nx,nx))
cd_s =np.zeros((n,2))
for day in days:
    file = 'trimdoz_'+day+'_01.nc'
    t_init  = 39600
    t_final = 97200
    pre_dt = datetime.datetime(int(file[8:12]),int(file[12:14]),int(file[14:16]),0,0)
    stdt = pre_dt+datetime.timedelta(seconds=t_init)
    endt = pre_dt+datetime.timedelta(seconds=t_final)
    Hg,Hv = read_sfc_data('sh',nt,stdt)
    H_s[i,:]=Hg[t,:]
    plt.figure()
    plt.hist(Hv[t,::2],bins=50)

# %%
H_s.shape
from numpy import random

x=random.randint(100, size=(5))


# %%
def checkerboard(sz=1,int_=5,plus=5):
    other=True
    a=np.zeros((20,20))
    for i in range(int(20/sz)):
        other= not other
        for j in range(int(20/sz)):
            if other:
                a[i*sz:i*sz+sz,j*sz:j*sz+sz]=random.randint(int_)+plus
                other=False
            else:
                a[i*sz:i*sz+sz,j*sz:j*sz+sz]=random.randint(int_)-plus
                other=True
    return a


# %%
H1=checkerboard(1,5,10)
H2=checkerboard(2,10,20)
H5=checkerboard(5,20,40)
H10=checkerboard(10,40,80)

# %%
plt.subplot(2,2,1)
plt.imshow(H1)
plt.title('H1:0-5')

plt.subplot(2,2,2)
plt.imshow(H2)
plt.title('H2:5-10')

plt.subplot(2,2,3)
plt.imshow(H5)
plt.title('H5:20-25')

plt.subplot(2,2,4)
plt.imshow(H10)
plt.title('H10:')

# %%
np.sqrt(2*25000**2)

# %%
cuts=np.linspace(0,.5)
h1s=np.zeros(cuts.shape)
h2s=np.zeros(cuts.shape)
h5s=np.zeros(cuts.shape)
h10s=np.zeros(cuts.shape)
for i in range(len(cuts)):
    cut=cuts[i]
    h1s[i]=el_het(-1,H1,cut)
    h2s[i]=el_het(-1,H2,cut)
    h5s[i]=el_het(-1,H5,cut)
    h10s[i]=el_het(-1,H10,cut)
plt.plot(cuts,h1s,alpha=.5)
plt.plot(cuts,h2s,alpha=.5)
plt.plot(cuts,h5s,alpha=.5)
plt.plot(cuts,h10s,alpha=.5)
plt.legend(['1','2','5','10'])

# %%
for i in range(len(days)):
    #fpzt=nc.Dataset('/home/tsw35/tyche/clubb/circ_tune2dlo2/'+days[i]+'/cr0.01k2fon1/k_2/agg_outzt.nc')
    plt.figure()
    plt.imshow(clsts[i,:],cmap='cool_r',origin='lower')
    plt.title(days[i]+'--'+str(cd_s[i,:]))
    plt.arrow(10,10,10*cd_s[i,0]/np.sum(cd_s[i,:]),10*cd_s[i,1]/np.sum(cd_s[i,:]),width=.5,length_includes_head=True)
    #plt.colorbar()
    plt.axis('off')
    break

# %%
cuts=np.array([.01,.05,.1,.2,.3])
lhets=np.zeros((len(days),len(cuts)))
for i in range(len(days)):
    for j in range(len(cuts)):
        lhets[i,j]=el_het(-1,clsts[i,:],cut=cuts[j])

# %%
binz

# %%
binz=np.linspace(2500,75000+2500,16)
for j in range(len(cuts)):
    plt.hist(lhets[:,j],histtype="stepfilled",alpha=.25,bins=binz)
plt.legend(cuts)

# %%
for i in np.where(lhets[:,1]<=25000)[0]:
    plt.figure()
    H=H_s[i,:]
    plt.imshow(H,vmin=np.percentile(H,2),vmax=np.percentile(H,98))

# %%

# %%
plt.hist(lhets[:,1],histtype="stepfilled",alpha=.5,bins=binz)
plt.hist(lhets[:,0],histtype="stepfilled",alpha=.5,bins=binz)

# %%
from scipy import interpolate


# %%
#### CALCULATE LENGTHSCALE OF HETEROGENEITY ####
def el_het(l_het,Hg_,cut=.25,samp=400):
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

        idx = np.random.choice(len(H_flat),size=samp,replace=True)
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
        #plt.plot(bins[0:-1],means)
        l_het_ = np.mean(bins[0:-1][means<=(cut*means[0])][0:1])
        #f= interpolate.interp1d(means,bins[0:-1]/2+bins[1:]/2)
        #l_het_ = f(cut*means[0])
        #plt.plot([bins[0],bins[-1]],[means[0]*cut,means[0]*cut],'k--',alpha=.5)
        #plt.title(l_het_)
        return l_het_

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
for day in days:
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    break

# %%
fpclst=nc.Dataset(dir2c+'k_2/clusters.nc','r')
fpzc  =nc.Dataset(dir2c+'k_2/c_1/output/arm_zt.nc','r') 
fpzh  =nc.Dataset(dir2c+'k_2/c_2/output/arm_zt.nc','r') 

# %%
for v in fpclst.variables:
    print(v)

# %%
m_dT=[]
m_stdT=[]
m_aT=[]
m_stdaT=[]
for day in days:
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    fpclst=nc.Dataset(dir2c+'k_2/clusters.nc','r')
    fpzc  =nc.Dataset(dir2c+'k_2/c_1/output/arm_zt.nc','r') 
    fpzh  =nc.Dataset(dir2c+'k_2/c_2/output/arm_zt.nc','r')
    fpm   =nc.Dataset(dir2c+'k_2/agg_outzm.nc','r')
    clst=fpclst['cluster'][0,:]
    tskin_clst=np.zeros((18,2))
    tskin=fpclst['tskin'][:]
    std_t=np.zeros((18,))
    del_t=np.zeros((18,))
    times_sfc=fpclst['t_model'][:]
    del_at=fpzh['thlm'][:,2,0,0]-fpzc['thlm'][:,2,0,0]
    std_at=2*np.sqrt(fpm['thlp2'][:,2,0,0])
    for t in range(18):
        std_t[t]=2*np.std(tskin[t,:])
        tskin_clst[t,0]=np.mean(tskin[t,:][clst==0])
        tskin_clst[t,1]=np.mean(tskin[t,:][clst==1])
        del_t[t]=tskin_clst[t,1]-tskin_clst[t,0]
        
    m_dT.append(del_t[int(len(std_t)/2)])
    m_stdT.append(std_t[int(len(std_t)/2)])
    m_aT.append(del_at[int(len(std_at)/2)])
    m_stdaT.append(std_at[int(len(std_at)/2)])
    
    plt.figure()
    plt.plot(np.linspace(5,22,18),std_t,'-o')
    plt.plot(np.linspace(5,22,18),del_t,'-o')
    plt.plot(np.linspace(5,22,len(del_at)),del_at)
    plt.plot(np.linspace(5,22,len(del_at)),std_at)


# %%
#m_dT=[]
#m_stdT=[]
#m_aT=[]
#m_stdaT=[]
binz=np.linspace(0,10,21)
x1=np.histogram(np.array(m_dT)*.8,binz)
plt.plot(binz[0:-1]+.25,x1[0],'-o')
x2=np.histogram(m_stdT,binz)
plt.plot(binz[0:-1]+.25,x2[0],'-o')
x3=np.histogram(m_aT,binz)
plt.plot(binz[0:-1]+.25,x3[0],'-o')
x4=np.histogram(m_stdaT,binz)
plt.plot(binz[0:-1]+.25,x4[0],'-o')
plt.legend(['T','stdT','aT','stdaT'])

# %%
fpsc=nc.Dataset(dir2c+'k_2/c_1/output/arm_sfc.nc','r')
fpsc=nc.Dataset(dir2c+'k_2/c_1/output/arm_sfc.nc','r') 

# %%
for v in fpsc.variables:
    print(v)

# %%

# %%
.75*1.5

# %%
c=np.array([.75,.25])
c=c/np.sqrt(np.sum(c**2))

# %%
c

# %%
c2=c*3
c3=c2-np.array([1.5,5])
c3[c3<0]=0

# %%
c3

# %%
