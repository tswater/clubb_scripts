import netCDF4 as nc
import shutil
import subprocess
import os
import datetime
import argparse
from mpi4py import MPI
import numpy as np

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ----------- #
# USER INPUTS #
# ----------- #
clubb_dir='/home/tsw35/tyche/clubb/'
search_dirs=['circ_tune']
sd=[]
for file in search_dirs:
    sd.append(clubb_dir+file)
# ---------------- #
# HELPER FUNCTIONS #
# ---------------- #

#### CALC_BLH ####
def calc_blh(ri_,alt_):
    blh=np.zeros((ri_.shape[0],))
    for t in range(ri_.shape[0]):
        
        try:
            wheres=np.where(ri_[t,:]>.25)[0]
            idx=wheres[0]
            if (idx<=3)&(len(wheres)>=2):
                idx=wheres[1]
            if (idx<=3)&(len(wheres)>=3):
                idx=wheres[2]
            if (idx<=3)&(len(wheres)>=4):
                idx=wheres[3]
            blh[t]=alt_[idx]
        except:
            pass
    return blh

#### FIND Ks ####
def find_ks(dir_):
    out=[]
    try:
        filelist=os.listdir(dir_)
    except:
        return []
    for file in filelist:
        if 'k_' in file:
            out.append(dir_+'/'+file)
        else:
            out.extend(find_ks(dir_+'/'+file))
    return out
# ----------- #
# CORE SCRIPT #
# ----------- #
k_dirs=[]
for file in sd:
    k_dirs.extend(find_ks(file))
i=0
for k_folder in k_dirs[rank::size]:
    print('RANK: '+str(rank)+'  '+str(i)+'/'+str(len(k_dirs)/size),flush=True)
    if (k_folder[-1]=='2'):
        fp1=nc.Dataset(k_folder+'/c_1/output/arm_zm.nc','r')
        fp2=nc.Dataset(k_folder+'/c_2/output/arm_zm.nc','r')
        blh_1=calc_blh(fp1['Richardson_num'][:,:,0,0],fp1['altitude'][:])
        blh_2=calc_blh(fp2['Richardson_num'][:,:,0,0],fp2['altitude'][:])
        fp1.close()
        fp2.close()
        fp1=nc.Dataset(k_folder+'/c_1/output/arm_sfc.nc','r+')
        fp2=nc.Dataset(k_folder+'/c_2/output/arm_sfc.nc','r+')
        dim=('time','altitude','latitude','longitude')
        try:
            fp1.createVariable('blh','d',dimensions=dim)
        except:
            pass
        try:
            fp2.createVariable('blh','d',dimensions=dim)
        except:
            pass
        fp1['blh'][:]=blh_1[:]
        fp2['blh'][:]=blh_2[:]
        fp1.close()
        fp2.close()
        fp3=nc.Dataset(k_folder+'/agg_outsfc.nc','r+')
        try:
            fp3.createVariable('blh','d',dimensions=dim)
        except:
            pass
        fp3['blh'][:]=(blh_1[:]+blh_2[:])/2
        fp3.close()
    elif (k_folder[-1]=='1'):
        fp1=nc.Dataset(k_folder+'/c_1/output/arm_zm.nc','r')
        blh_1=calc_blh(fp1['Richardson_num'][:,:,0,0],fp1['altitude'][:])
        fp1.close()
        fp1=nc.Dataset(k_folder+'/c_1/output/arm_sfc.nc','r+')
        try:
            fp1.createVariable('blh','d',dimensions=dim)
        except:
            pass
        fp1['blh'][:]=blh_1[:]
        fp1.close()
    i=i+1








