import netCDF4 as nc
import shutil
import subprocess
import os
import datetime
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
cpl_dir = '/home/tsw35/tyche/clubb/sgp_cpl/'
sgp_dir = '/stor/soteria/hydro/shared/lasso_for_tyler/'

def strptime(dtstr):
    return datetime.datetime(int(dtstr[0:4]),int(dtstr[5:7]),int(dtstr[8:10]),\
           int(dtstr[11:13]),int(dtstr[14:16]))



# --------- #
# CORE CODE #
# --------- #
# create basic directory
try:
    os.mkdir(cpl_dir)
except:
    pass

# create a list of all valid folders for sgp runs
allfiles = []
for file in os.listdir('../'):
    if ('sgp' not in file) or ('cpl' in file):
        continue
    else:
        allfiles.append(file)
print(len(allfiles))
# iterate through them
for folder in os.listdir(sgp_dir)[rank::size]:
    start_time = datetime.datetime(2022,9,9,9,9)
    end_time   = datetime.datetime(2000,2,2,2,2)
    for file in os.listdir(sgp_dir+'/'+folder):
        if len(file)<25:
            continue
        time = strptime(file[-16:])
        if time < start_time:
            start_time = time
        elif time > end_time:
            end_time = time
    file="sgp_"+str(start_time.year)+str(start_time.month)+str(start_time.day)


    # find the date and time of the run
    fp = open('/home/tsw35/tyche/clubb/'+file+'/k_1/c_1/input/arm_model.in','r')
    for line in fp.readlines():
        if line[0:3]=='day':
            lp = line.split(' ')
            lp = [x for x in lp if x !='']
            day = int(lp[2])
        elif line[0:5]=='month':
            lp = line.split(' ')
            lp = [x for x in lp if x !='']
            month = int(lp[2])
        elif line[0:4]=='year':
            lp = line.split(' ')
            lp = [x for x in lp if x !='']
            year = int(lp[2])
        elif line[0:12]=='time_initial':
            lp = line.split(' ')
            lp = [x for x in lp if x !='']
            t_init = float(lp[2])
        elif line[0:10]=='time_final':
            lp = line.split(' ')
            lp = [x for x in lp if x !='']
            t_final = float(lp[2])
    fp.close()
    pre_dt = datetime.datetime(year,month,day,0,0)
    stdt = pre_dt+datetime.timedelta(seconds=t_init)
    endt = pre_dt+datetime.timedelta(seconds=t_final)
    
    # convert to iso strings
    st_str = stdt.strftime('%Y-%m-%dT%H:%M:%S')
    en_str = endt.strftime('%Y-%m-%dT%H:%M:%S')
    
    # run the simple script
    #cmd1 = 'python cpl_lrg_frc.py -k '+str(k)+' -i '+cpl_dir[3:]+file+' -s '+st_str+' -e '+en_str+' -c '+sgp_dir+folder+'/ -n '+str(nx)+' -d '+str(dx)
    cmd1 = 'python cpl_ant.py -l 40000 -d '+str(dx)+' -i '+cpl_dir[-8:]+file+' -s '+st_str+' -e '+en_str+' -c '+sgp_dir+folder+'/ -n '+str(nx)+' -k '+str(k)
    print('#######\n'+cmd1+'\n'+'rank: '+str(rank)+'\n#######')
    cmd2 = 'python cpl_agg.py -i '+cpl_dir+file+'/'
    subprocess.run(cmd1,shell=True)
    subprocess.run(cmd2,shell=True)
    print('#######\n'+cmd2+'\n'+'rank: '+str(rank)+'\n#######')


