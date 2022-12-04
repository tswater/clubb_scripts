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
k = 2
dx = 250
nx = 401
delta_t = 30
lhet = -2
flux_on = 1
wind_cr = 1
sgp_dir  = '/stor/soteria/hydro/shared/lasso_for_tyler/'
day_dir  = '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
sgp_dir2 = '/stor/soteria/hydro/shared/data/tylersclutter/surfaces201589/'
agg_only = False

# Lits
#cratios=[10,15,20,25,30]
cratios=[20]
#cts    =[.0001,.00015,.0002,.00025,.0003,.00035,.0004]
cutoffs=[.8]
#cutoffs=[.5,.6,.7,.8,.9]
cts     =[.00001,.00002,.00003,.00004,.00005,.00006,.00007,.00008,.00009,\
          .00041,.00042,.00043,.00044,.00045,.00046,.00047,.00048,.00049,.0005]
          #np.linspace(.0001,.0004,31)

#days = ['20180805','20150829','20160719','20170716','20170705','20180916','20190804','20160619','20180710','20170923']
days=[]
for file in os.listdir(day_dir):
    if '_01' in file:
        days.append(file[8:16])
    else:
        continue

runlist=[]
for day in days:
    # make folder for the day
    try:
        os.mkdir(clubb_dir+'circ_tune/'+day)
    except:
        pass
    for cr in cratios:
        for ct in cts:
            for cf in cutoffs:
                runlist.append([day,cr,ct,cf])
if rank==0:
    print(len(runlist))

for run in runlist[rank::size]:
    file = 'trimdoz_'+run[0]+'_01.nc'
    datest = file[8:12]+'-'+file[12:14]+'-'+file[14:16]+'-'+'15'
    try:
        sfc_dir=sgp_dir+'sgp_'+file[8:16]+'_00/'
        sfc_file=sfc_dir+'jsssh_bdy_02_'+datest+'-00'
        fp=open(sfc_file,'r')
        fp.close()
    except:
        sfc_dir=sgp_dir2
        sfc_file=sfc_dir+'jsssh_bdy_02_'+datest+'-00'
        fp=open(sfc_file,'r')
        fp.close()

    t_init  = 43200
    t_final = 100800
    pre_dt = datetime.datetime(int(file[8:12]),int(file[12:14]),int(file[14:16]),0,0)
    stdt = pre_dt+datetime.timedelta(seconds=t_init)
    endt = pre_dt+datetime.timedelta(seconds=t_final)

    # convert to iso strings
    st_str = stdt.strftime('%Y-%m-%dT%H:%M:%S')
    en_str = endt.strftime('%Y-%m-%dT%H:%M:%S')
    
    runstr='cr'+str(run[1])+'ct'+str(run[2])+'cf'+str(run[3])

    cmd1 = 'python cpl_ant.py -l '+str(lhet)+' -d '+str(dx)+\
           ' -i '+'circ_tune/'+run[0]+'/'+runstr+' -s '+st_str+' -e '+en_str+\
           ' -c '+sfc_dir+'/ -n '+str(nx)+' -k '+str(k)+\
           ' -f '+str(flux_on)+' -w '+str(wind_cr)+\
           ' -y '+str(delta_t)+' -r '+str(run[1])+\
           ' -t '+str(run[2])+' -u '+str(run[3])

    aggdir=clubb_dir+'circ_tune/'+run[0]+'/'+runstr+'/'
    if not agg_only:
        print('#######\n'+cmd1+'\n'+'rank: '+str(rank)+'\n#######')
        cmd2 = 'python cpl_agg.py -i '+aggdir
        subprocess.run(cmd1,shell=True)
    cmd2 = 'python cpl_agg.py -i '+aggdir
    print('#######\n'+cmd2+'\n'+'rank: '+str(rank)+'\n#######')
    subprocess.run(cmd2,shell=True)


