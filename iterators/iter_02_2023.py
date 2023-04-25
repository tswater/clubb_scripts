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
delta_t = 30
lhet = -1
flux_on = 1
wind_cr = 1
sgp_dir  = '/stor/soteria/hydro/shared/lasso_for_tyler/'
day_dir  = '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
sgp_dir2 = '/stor/soteria/hydro/shared/data/tylersclutter/surfaces201589/'
agg_only = False
k1only   = True

os.chdir('/home/tsw35/soteria/clubb/clubb_scripts')

mainfolder= 'tune_feb23/'

# Lits
curs = [.5,1,2,3.5,5,7,10]
meth = ['denv_ds','denv_s','denv_a']
cc1s = [.7,.8,.85,.9,.95,1]
cc2s = [1]

#days = ['20160716','20180619','20160719','20150829','20170716','20170705','20180916','20190707','20190714','20170923','20160625','20170627','20170728','20180707','20180916']
days = ['20160716','20180619','20160719','20150829','20170716','20190707','20190714','20160625','20170627','20180707']
days.sort()

runlist=[]
for day in days:
    # make folder for the day
    try:
        os.mkdir(clubb_dir+mainfolder+day)
    except Exception as e:
        print(e)
        pass
    for cr in curs:
        for cc in cc1s:
            for inc_ in cc2s:
                for mt in meth:
                    if k1only:
                        continue
                    runlist.append([day,cr,2,1,cc,inc_,mt])
    #runlist.append([day,0,2,0,0,0,'h_h'])
    runlist.append([day,0,1,0,0,0,'h_h'])
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

    t_init  = 36000
    t_final = 97200
    pre_dt = datetime.datetime(int(file[8:12]),int(file[12:14]),int(file[14:16]),0,0)
    stdt = pre_dt+datetime.timedelta(seconds=t_init)
    endt = pre_dt+datetime.timedelta(seconds=t_final)

    # convert to iso strings
    st_str = stdt.strftime('%Y-%m-%dT%H:%M:%S')
    en_str = endt.strftime('%Y-%m-%dT%H:%M:%S')
    
    runstr='feb23_cr'+str(run[1])+'k'+str(run[2])+'cc'+str(run[4])+'mm'+str(run[6])

    cmd1 = 'python cpl_ant.py -l '+str(lhet)+\
           ' -i '+mainfolder+run[0]+'/'+runstr+' -s '+st_str+' -e '+en_str+\
           ' -k '+str(run[2])+ ' -c2 '+str(run[5])+' -c1 '+str(run[4])+\
           ' -w '+str(wind_cr)+' -cm '+str(run[6])+\
           ' -y '+str(delta_t)+' -r '+str(run[1])+' -f '+str(run[3])
    aggdir=clubb_dir+mainfolder+run[0]+'/'+runstr+'/'
    if not agg_only:
        print('#######\n'+cmd1+'\n'+'rank: '+str(rank)+'\n#######')
        cmd2 = 'python cpl_agg.py -i '+aggdir
        subprocess.run(cmd1,shell=True)
    cmd2 = 'python cpl_agg.py -i '+aggdir
    print('#######\n'+cmd2+'\n'+'rank: '+str(rank)+'\n#######')
    subprocess.run(cmd2,shell=True)

