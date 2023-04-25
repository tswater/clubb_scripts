import numpy as np
import netCDF4 as nc
from sklearn.linear_model import LinearRegression
import subprocess
import datetime
import os
import argparse
from mpi4py import MPI

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

subdir='apr25_dz60'
les_dir='/home/tsw35/soteria/clubb/data/les_param2/'
dt    = 6
dz    = 60
prs = argparse.ArgumentParser(description='Short sample app')
prs.add_argument('-d', action='store', dest='dr', default=subdir)
prs.add_argument('-z', action='store', dest='dz',type=int, default=dz)
prs.add_argument('-t', action='store', dest='dt',type=float, default=dt)

args = prs.parse_args()
subdir=args.dr
dz=args.dz
dt=args.dt

cpl_dir  = '/home/tsw35/tyche/clubb/'+subdir

os.chdir('/home/tsw35/soteria/clubb/clubb_scripts')

filelist=os.listdir(les_dir)
filelist.sort()

bad_days=['20180811.nc','20160819.nc','20170705.nc','20190804.nc','20190805.nc']

cmds=[]
agg_cmds=[]

try:
    os.mkdir(cpl_dir+'_cpl')
except:
    pass
try:
    os.mkdir(cpl_dir+'_2c')
except:
    pass
try:
    os.mkdir(cpl_dir+'_1c')
except:
    pass

filelist=['20150801.nc','20160610.nc','20160625.nc','20160716.nc','20170716.nc','20170717.nc','20170802.nc','20180707.nc','20180709.nc','20190707.nc']

for file in filelist:
    if file in bad_days:
        continue
    fp=nc.Dataset(les_dir+file,'r')

    t1=3
    t2=12
    ht=16 #first 500m is 17
    ht0=3
    u90=fp['u90'][t1:t2,ht0:ht]
    u10=fp['u10'][t1:t2,ht0:ht]
    v90=fp['v90'][t1:t2,ht0:ht]
    v10=fp['v10'][t1:t2,ht0:ht]
    lhet=fp['lhet'][0]
    dvpt=fp['vpt'][t1:t2,1,ht0:ht]-fp['vpt'][t1:t2,0,ht0:ht]

    msk=((u90/np.abs(u90))==(u10/np.abs(u10)))|(dvpt<0)
    u90[msk]=0
    u10[msk]=0
    m2=(np.abs(u10)>np.abs(u90))
    u90[m2]=u10[m2]
    u90=np.abs(u90)
    
    msk=((v90/np.abs(v90))==(v10/np.abs(v10)))|(dvpt<0)
    v90[msk]=0
    v10[msk]=0
    m2=(np.abs(v10)>np.abs(v90))
    v90[m2]=v10[m2]
    v90=np.abs(v90)
    
    u2=np.sqrt(v90**2+u90**2)
    dvpt=dvpt[u2>0]
    u2  =u2[u2>0]
    
    if len(u2)<10:
        cf=1.3
        sc=-1
    else:
        um=dvpt/300*(lhet**(1/2)*9.8**(1/2))
        
        model=LinearRegression(fit_intercept=False)
        model.fit(um.reshape(-1, 1),u2.reshape(-1, 1))
        cf=model.coef_[0][0]
        sc=model.score(um.reshape(-1, 1),u2.reshape(-1, 1))
    if rank==0:
        print(file+': '+str(cf)+' : '+str(sc)) 

    t_init  = 36000+3600*2
    t_final = 97200+3600
    pre_dt = datetime.datetime(int(file[0:4]),int(file[4:6]),int(file[6:8]),0,0)
    stdt = pre_dt+datetime.timedelta(seconds=t_init)
    endt = pre_dt+datetime.timedelta(seconds=t_final)

    # convert to iso strings
    st_str = stdt.strftime('%Y-%m-%dT%H:%M:%S')
    en_str = endt.strftime('%Y-%m-%dT%H:%M:%S')
    
    cmd1 = 'python cpl_ant.py -i '+subdir+'_cpl/sgp_'+file[0:8]+' -s '+st_str+' -e '+en_str+\
           ' -y '+str(dt)+' -k 2 -f 1 -r '+str(cf)+' -z '+str(dz)
    cmd1a = 'python cpl_agg.py -i '+cpl_dir+'_cpl/sgp_'+file[0:8]+'/'
    
    cmd2 = 'python cpl_ant.py -i '+subdir+'_2c/sgp_'+file[0:8]+' -s '+st_str+' -e '+en_str+\
           ' -y '+str(dt)+' -k 2 -f 0'+' -z '+str(dz)
    cmd2a = 'python cpl_agg.py -i '+cpl_dir+'_2c/sgp_'+file[0:8]+'/'

    cmd3 = 'python cpl_ant.py -i '+subdir+'_1c/sgp_'+file[0:8]+' -s '+st_str+' -e '+en_str+\
           ' -y '+str(dt)+' -k 1 -f 0'+' -z '+str(dz)
    cmd3a = 'python cpl_agg.py -i '+cpl_dir+'_1c/sgp_'+file[0:8]+'/'
    
    cmds.append(cmd1)
    cmds.append(cmd2)
    cmds.append(cmd3)
    agg_cmds.append(cmd1a)
    agg_cmds.append(cmd2a)
    agg_cmds.append(cmd3a)

for i in list(range(len(cmds)))[rank::size]:
    print('#######\n'+cmds[i]+'\n'+'rank: '+str(rank)+'\n#######')
    subprocess.run(cmds[i],shell=True)
    print('#######\n'+agg_cmds[i]+'\n'+'rank: '+str(rank)+'\n#######')
    subprocess.run(agg_cmds[i],shell=True)

