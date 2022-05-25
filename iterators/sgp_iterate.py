import netCDF4 as nc
import shutil
import subprocess
import os
import datetime
import argparse
from mpi4py import MPI

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ----------- #
# USER INPUTS #
# ----------- #
k = 2
dx = 250
nx = 401
delta_t = .5 
lhet = -2
flux_on = 1
wind_cr = 1
ct = .00025
cpl_subdir = 'sgp_deleteme/'
sgp_dir  = '/stor/soteria/hydro/shared/lasso_for_tyler/'
day_dir  = '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
sgp_dir2 = '/stor/soteria/hydro/shared/data/tylersclutter/surfaces201589/'
agg_only = False # if True, will only aggregate not rerun model

#### ARG PARSER ####
prs = argparse.ArgumentParser(description='Short sample app')
prs.add_argument('-k', action='store', dest='k',type=int, default=k)
prs.add_argument('-l', action='store', dest='lh',type=int, default=lhet)
prs.add_argument('-w', action='store', dest='wc',type=int, default=wind_cr)
prs.add_argument('-f', action='store', dest='fn',type=int, default=flux_on)
prs.add_argument('-d', action='store', dest='dr', default=cpl_subdir)
prs.add_argument('-y', action='store', dest='dt',type=float, default=delta_t)
prs.add_argument('-a', action='store', dest='at', default='')
prs.add_argument('-t', action='store', dest='ct', type=float,default=ct)

args = prs.parse_args()
lhet    = args.lh
wind_cr = args.wc
flux_on = args.fn
k       = args.k
cpl_subdir = args.dr
delta_t = args.dt
atm_dir = args.at
ct      = args.ct

cpl_dir  = '/home/tsw35/tyche/clubb/'+cpl_subdir


def strptime(dtstr):
    return datetime.datetime(int(dtstr[0:4]),int(dtstr[5:7]),int(dtstr[8:10]),\
           int(dtstr[11:13]),int(dtstr[14:16]))

print('doing it',flush=True)


# --------- #
# CORE CODE #
# --------- #
# create basic directory
try:
    os.mkdir(cpl_dir)
except:
    pass

day_list=[]

for file in os.listdir(day_dir):
    if file[18]=='1':
        day_list.append(file)


for file in day_list[rank::size]:
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

    
    # run the simple script
    #cmd1 = 'python cpl_lrg_frc.py -k '+str(k)+' -i '+cpl_dir[3:]+file+' -s '+st_str+' -e '+en_str+' -c '+sgp_dir+folder+'/ -n '+str(nx)+' -d '+str(dx)
    cmd1 = 'python cpl_ant.py -l '+str(lhet)+' -d '+str(dx)+\
           ' -i '+cpl_subdir+'sgp_'+file[8:16]+' -s '+st_str+' -e '+en_str+\
           ' -c '+sfc_dir+'/ -n '+str(nx)+' -k '+str(k)+\
           ' -f '+str(flux_on)+' -w '+str(wind_cr)+\
           ' -y '+str(delta_t)+' -a '+str(atm_dir)+' -t'+str(ct)
    cmd2 = 'python cpl_agg.py -i '+cpl_dir+'sgp_'+file[8:16]+'/'
    if agg_only:
        print('#######\n'+cmd2+'\n'+'rank: '+str(rank)+'\n#######')
        subprocess.run(cmd2,shell=True)
    else:
        print('#######\n'+cmd1+'\n'+'rank: '+str(rank)+'\n#######')
        subprocess.run(cmd1,shell=True)
        print('#######\n'+cmd2+'\n'+'rank: '+str(rank)+'\n#######')
        subprocess.run(cmd2,shell=True)


