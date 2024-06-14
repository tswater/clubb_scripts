import numpy as np
import netCDF4 as nc
import scipy as sci
import scipy.ndimage
import os
#import pickle
import datetime
import rasterio

from mpi4py import MPI

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


#######################
#######################
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

############################################
############################################


odir='/home/tsw35/soteria/clubb/data/les_param2/'
lesdir='/home/tsw35/xTyc_shared/clasp/fr2/'
sigma=20
dirlist=os.listdir(lesdir)
dirlist.sort()
dthv={}
vpt={}
dtsk={}
lhet={}
dirlist0=[]
dirlist1=[]

u90={}
u10={}
v90={}
v10={}

for day in dirlist:
    if '_01' in day:
        continue
    dirlist0.append(day)
for day in dirlist0[rank::size]:
    dirlist1.append(day)

print(len(dirlist1))
for day in dirlist1:
    print(day,end=':')
    dst=day[4:12]
    dirlist2=os.listdir(lesdir+day)
    dirlist2.sort()
    nt=len(dirlist2)
    dthv[dst]=np.zeros((nt,2))
    vpt[dst]=np.zeros((nt,2,226))
    dtsk[dst]=np.zeros((nt,2))
    lhet[dst]=0
    u90[dst]=np.zeros((nt,226))
    u10[dst]=np.zeros((nt,226))
    v90[dst]=np.zeros((nt,226))
    v10[dst]=np.zeros((nt,226))



    #### Calculate Lengthscale of Heterogeneity
    hggt=5
    stdt=datetime.datetime(int(day[4:8]),int(day[8:10]),int(day[10:12]),12)
    try:
        Hg,Hv = read_sfc_data('sh',nt,stdt)
    except Exception as e:
        print(day)
        print(e)
        continue
    Hgg = Hg[hggt,:,:]
    lhet[dst] = estimate_l_het(-1,Hgg)
    
    for t in range(nt):
        print(t,end=',',flush=True)
        try:
            fname=lesdir+day+'/'+dirlist2[t]
            fpht=nc.Dataset(fname,'r')
            day2=day.replace('_00','_01')
            fname2=lesdir+day2+'/'+dirlist2[t]
            fphm=nc.Dataset(fname2,'r')
        except Exception as e:
            print(day)
            print(e)
            continue
        uht=fpht['AVV_U'][0,:,:,:]
        vht=fpht['AVV_V'][0,:,:,:]
        uhm=fphm['AVV_U'][0,:,:,:]
        vhm=fphm['AVV_V'][0,:,:,:]
        #H = fpht['AVS_SH'][0,:,:]
        th = fpht['AVV_THV'][0,:,:,:]
        tsk = fpht['AVS_TSK'][0,:,:]
        alt=fpht['AVP_Z'][0,:]

        dtsk[dst][t]=2*np.std(tsk)
        
        thvg = np.zeros(th.shape)

        for z in range(uht.shape[0]):
            
            uhtn=sci.ndimage.filters.gaussian_filter(uht[z,:],sigma)
            vhtn=sci.ndimage.filters.gaussian_filter(vht[z,:],sigma)

            uhmn=sci.ndimage.filters.gaussian_filter(uhm[z,:],sigma)
            vhmn=sci.ndimage.filters.gaussian_filter(vhm[z,:],sigma)
            
            u10[dst][t,z]=np.percentile(uhtn,10)-np.percentile(uhmn-np.mean(uhmn),10)
            u90[dst][t,z]=np.percentile(uhtn,90)-np.percentile(uhmn-np.mean(uhmn),90)

            v10[dst][t,z]=np.percentile(vhtn,10)-np.percentile(vhmn-np.mean(vhmn),10)
            v90[dst][t,z]=np.percentile(vhtn,90)-np.percentile(vhmn-np.mean(vhmn),90)

            thvg[z,:,:]=sci.ndimage.filters.gaussian_filter(th[z,:,:],10)
            
        dthv[dst][t]=2*np.std(thvg[10,:,:])
        
        clst=np.zeros(thvg[0,:,:].shape)
        clst[thvg[10,:,:]<np.percentile(thvg[10,:,:],50)]=1
        clst[thvg[10,:,:]>np.percentile(thvg[10,:,:],50)]=2
        if np.all(clst==0) or np.all(clst<=1):
            vpt[dst][t,0,:]=np.mean(thvg[0:226,:,:],axis=(1,2))
            vpt[dst][t,1,:]=np.mean(thvg[0:226,:,:],axis=(1,2))
        else:
            vpt[dst][t,0,:]=np.mean(thvg[0:226,clst==1],axis=1)
            vpt[dst][t,1,:]=np.mean(thvg[0:226,clst==2],axis=1)

    print('...',end='')
    fpo=nc.Dataset(odir+dst+'.nc','w')
    fpo.createDimension('z',size=226)
    fpo.createVariable('z','d',dimensions=('z'))
    fpo['z'][:]=alt[:]
    fpo.createDimension('clst',size=2)
    fpo.createDimension('time',size=nt)
    fpo.createVariable('vpt','d',dimensions=('time','clst','z'))
    fpo.createVariable('dthv','d',dimensions=('time','clst'))
    fpo.createVariable('dtsk','d',dimensions=('time','clst'))
    fpo.createVariable('lhet','d',dimensions=('time'))
    fpo.createVariable('u10','d',dimensions=('time','z'))
    fpo.createVariable('u90','d',dimensions=('time','z'))
    fpo.createVariable('v10','d',dimensions=('time','z'))
    fpo.createVariable('v90','d',dimensions=('time','z'))

    fpo['vpt'][:]=vpt[dst][:]
    fpo['dtsk'][:]=dtsk[dst][:]
    fpo['dthv'][:]=dthv[dst][:]
    fpo['lhet'][:]=lhet[dst]
    fpo['u10'][:]=u10[dst][:]
    fpo['u90'][:]=u90[dst][:]
    fpo['v10'][:]=v10[dst][:]
    fpo['v90'][:]=v90[dst][:]
    fpo.close()
    print()

#fp=open('lesmar3.p','wb')
#pickle.dump(uu,fp)
#pickle.dump(dthv,fp)
#pickle.dump(vpt,fp)
#pickle.dump(dtsk,fp)
#pickle.dump(lhet,fp)

#fp.close()

